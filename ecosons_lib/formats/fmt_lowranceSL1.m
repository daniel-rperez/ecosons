%%%fmt_lowranceSL1.m%%%
%Reads a single-beam Lowrance sonar file into memory
% [P HS PS]=fmt_lowranceSL1(fname)
% P{}: channel cells with ping matrices: rows: ping no., columns: echo sample
% HS: channel cells with transducer headers
% PS: GPS + time data (time=-1: no data)
% fname: RAW file filename
function [P,HS,PS]=fmt_lowranceSL1(fname)
 %l2=log2(10); %%%in case file is codified as binary logarithms
 %fname = strrep(fname, "\\", "/");

 %Initialize returned values
 P=[];
 HS=struct;
 PS=struct;
  lPS=gpsLowrance(NaN,NaN,NaN); %%%function to read and decodified coordinates information

 %file time
 s=stat(fname);
 ftime=localtime(s.mtime); %%%keep and convert local time of adquisition

 %open file
 sonarFile=fopen(fname, 'rb');

 %compute file length
 fseek(sonarFile,0,"eof");
 sonarFileSize=ftell(sonarFile);
 fseek(sonarFile,0,"bof");

 %channel number?
 nchan=fread(sonarFile, 1, "int32"); %%%ask number of available channels

 %block size
 blksz=fread(sonarFile, 1, "uint8")+...%%%read and conversion information to integer numbers
       256*fread(sonarFile, 1, "uint8")+...
       65536*fread(sonarFile, 1, "uint8");

 %compression scheme (only schem=0 is currently read)
 schem=fread(sonarFile, 1, "uint8");

 %number of pings in the file
 npings=floor(sonarFileSize/blksz);

 %reserve space
 P=zeros(npings, blksz-3*4); %%%new matriz to keep bin information with number of columns equal to blksz size without three values previous to bins information (cmask, lowerLimit and depth)

 %%%define and initiate variables
 X=NaN;
 Y=NaN;
 T=0;
 temperature=0;
 soundSpeed=1500;
 T0=(ftime.hour+ftime.min/60+ftime.sec/3600+ftime.usec/3600e+6); %%%calculate file recorded decimal hour (finished it)

 %read pings
 for n=1:npings

  %field bitmask
  vmask=fread(sonarFile, 1, "uint32");

  lowerLimit=0.3048*fread(sonarFile, 1, "float32"); %%%read lowerLimit of each ping turn from feet to meters
  depth=0.3048*fread(sonarFile, 1, "float32"); %%%read ping depth turn from feet to meters

  hdrsz=3*4;

    # {'name': 'surface_depth', 'type': 'f'},
    # {'name': 'top_of_bottom', 'type': 'f'},
    # {'name': 'temperature2', 'type': 'f'},
    # {'name': 'temperature3', 'type': 'f'},
    # {'name': 'time1', 'type': 'I'},
    # {'name': 'gps_speed', 'type': 'f'},
    # {'name': 'heading', 'type': 'f'},
    # {'name': 'altitude', 'type': 'f'},
    # {'name': 'packetsize', 'type': 'H'},

  %assumed to be zero
  if( bitand(vmask, 0x00080000 ) )
   hdrsz=hdrsz+4;
   upperLimit=0.3048*fread(sonarFile,1,"float32");
  else
   upperLimit=0;
  endif

  %not currently used
  if( bitand(vmask, 0x00100000 ) )
   hdrsz=hdrsz+4;
   temperature=0.3048*fread(sonarFile,1,"float32");
  else
   temperature=0;
  endif

  %not currently used
  if( bitand(vmask, 0x00800000 ) )
   hdrsz=hdrsz+4;
   waterSpeed=0.514444*fread(sonarFile,1,"float32");
  else
   waterSpeed=NaN;
  endif

  %GPS coordinates
  if( bitand(vmask, 0x01000000 ) )
   hdrsz=hdrsz+4;
   Y=fread(sonarFile,1,"uint32");
   hdrsz=hdrsz+4;
   X=fread(sonarFile,1,"uint32");
  endif

  %not currently used
  if( bitand(vmask, 0x04000000 ) )
   hdrsz=hdrsz+4;
   surfaceDepth=0.3048*fread(sonarFile, 1, "float32");
  endif
  if( bitand(vmask, 0x08000000 ) )
   hdrsz=hdrsz+4;
   topOfBottomDepth=0.3048*fread(sonarFile, 1, "float32");
  else
   topOfBottomDepth=NaN;
  endif

  %not used
  if( bitand(vmask, 0x00200000 ) )
   hdrsz=hdrsz+4;
   temperature2=fread(sonarFile,1,"float32");
  endif
  if( bitand(vmask, 0x00400000 ) )
   hdrsz=hdrsz+4;
   temperature3=fread(sonarFile,1,"float32");
  endif

  %time in milliseconds
  if( bitand(vmask, 0x20000000 ) )
   hdrsz=hdrsz+4;
   T=fread(sonarFile, 1, "uint32");
  endif
  
  %not used
  if( bitand( vmask, 0x40000000 ) )
   hdrsz=hdrsz+4;
   gps_speed=fread(sonarFile,1,"float32");
   hdrsz=hdrsz+4;
   heading=fread(sonarFile,1,"float32");
  endif
  if( bitand( vmask, 0x00040000 ) )
   hdrsz=hdrsz+4;
   altitude=fread(sonarFile,1,"float32");
  endif
  if( bitand( vmask, 0x10000000 ) )
   hdrsz=hdrsz+4;
   packetsize=fread(sonarFile,1,"uint16");
  endif

  %load ping
  if( schem == 0 )
   %P(n, 1:(blksz-hdrsz))=l2*fread(sonarFile, blksz-hdrsz, "uint8"); %in case 2-log base codification (not proved)
   P(n, 1:(blksz-hdrsz))=fread(sonarFile, blksz-hdrsz, "uint8"); %%%read bin information, from current point to final echo information (block length minus head length)
  else
   %P(n, 1:(blksz-hdrsz))=fread(sonarFile, blksz-hdrsz, "uint8"); %%%no schem!=0 has been tested
   fseek(f, blksz-hdrsz, "cof");
   P(n, 1:(blksz-hdrsz))=nan(1,blksz-hdrsz);
  endif

  %assume data is right
  %%%calculate bin length (in time)
  binLength=2*lowerLimit/(soundSpeed*(blksz-hdrsz)); %%%bin length in meters

  %create structure (fields inherited from Simrad's format; not all present in Lowrance SLG)
  samp=struct;
   samp.channel=nchan;
   samp.frequency=200e+3; %%%usual carrier frequency; not included in the file
   samp.transmitPower=1000; %%%place holder for actual emission power (unknown)
   samp.sampleInterval=binLength;
   samp.pulseLength=4*binLength;
   samp.soundVelocity=soundSpeed;
   samp.absorptionCoefficient=0.0;
   samp.temperature=temperature;
   samp.count=n;
   %store some other data (mainly for debugging purposes)
   samp._depth=depth; %%%informative or computed by the echosounder variables; dispensable ones
   samp._tobDepth=topOfBottomDepth;
   samp._schem=schem;

  HS(n)=samp;
  PS(n)=gpsLowrance(X,Y,T);

 endfor

 %set time in hours
 T0=T0-PS(end).time/3600;
 for n=1:npings
  PS(n).time=PS(n).time/3600+T0;
 endfor

 %return data as cell
 P={ P };
 HS={ HS };

endfunction

%%%GPS Coordinate conversion
function lPS=gpsLowrance(X,Y,T)
 Sm=6356752.3142; %terrestrial radius
 dg=57.2957795132; %180/pi;

 lPS._X=X;
 lPS._Y=Y;

 lPS.time=-1.0;
 lPS.latitude=0.0;
 lPS.longitude=0.0;

 X=bitand(0x00FFFFFF,X);
 if( bitand(0x00800000,X) )
  X=double(X)-(256^3-1);
 endif

 Y=bitand(0x00FFFFFF,Y);
 if( bitand(0x00800000,Y) )
  Y=double(Y)-(256^3-1);
 endif

 lPS.latitude=dg*(2*atan(exp(Y/Sm))-pi/2);
 lPS.longitude=dg*X/Sm;
 lPS.time=T/1000;
endfunction

