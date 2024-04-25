%Reads a single-beam Lowrance sonar file (v.2) into memory
% [P Q G]=fmt_lowranceSL2(fname)
% P{}: channel cells with ping matrices: rows: ping no., columns: echo sample
% Q{}: channel cells with transducer headers
% G: GPS + time data (time=-1: no data)
% fname: SL2 file filename
function [P,Q,G]=fmt_lowranceSL2(fname)

  Sm=6356752.3142; %radio terrestre
  dg=57.2957795132; %180/pi;
  temperature=0;
  soundSpeed=1500;

	f=fopen(fname, "rb");

	%channel codes
	nchannels=[];
	%number of pings per channel
	npings=[];
	%maximum echo length per channel
	lpings=[];
	%total number of pings
	ncount=0;

	%hd0
	hd0=fread(f,8,'uint8');
	%packet size (datagram+echo)
	ll=fread(f,1,'uint32');

	while( ~feof(f) )
		%datagram header
		hd=fread(f, 140, "uint8");

		nc=hd(29); %channel number
		if( ~any( nc == nchannels ) )
			nchannels=[nchannels, nc];
			npings=[npings, 0];
			lpings=[lpings, 0];
		endif

    %ping number (in channel?); keep maximum
    n=typecast(uint8(hd(33:36)'), "uint32")+1;
    npings(nchannels==nc)=max(n, npings(nc == nchannels));

		%echo length; keep maximum
		l=hd(31)+256*hd(32);
    lpings(nc == nchannels)=max(l, lpings(nc == nchannels));

		%skip echo
		fseek(f, l, "cof");

		%check consistency
		ll=fread(f,1,'uint32');
		if( ll+4 ~= ftell(f) )
			fclose(f);
      error('Error in file');
			return
		endif

		%count datagrams
		ncount=ncount+1;
	endwhile

  %pings and headers lists
	P=[];
	Q=struct;
  G=struct;

  %GPS data
	X=nan(1,ncount);
	Y=nan(1,ncount);
	T0=-1;
	T=nan(1,ncount);
	D=nan(1,ncount);

	for nc=nchannels(1)

		%back to the beginning
		fseek(f, 0, "bof");

		%hd0
		hd0=fread(f,8,'uint8');
		ll=fread(f,1,'uint32');

		nping=0;
		for n=1:ncount
			%datagram header
			hd=fread(f, 140, "uint8");
      l=hd(31)+256*hd(32);

      if( T0 < 0 )
	  		T0=double( typecast(uint8(hd(57:60)'), "uint32") );
		  endif

 			if( nc == hd(29) )
        nping=nping+1;

  			%GPS
        Gn=struct;
        xy=typecast(uint8(hd(105:112)'), "int32");
        Gn.longitude=dg*double(xy(1))/Sm;
        Gn.latitude=dg*(2*atan(exp(double(xy(2))/Sm))-pi/2);
	  		tt=typecast(uint8(hd(137:140)'), "int32"); %??? datos de tempo GPS en ms dende o comezo do arquivo?
        dte=localtime(T0+double( tt )/1000);
        Gn.time=dte.hour+(dte.min+(dte.sec+dte.usec/1000000)/60)/60;
        Gn._year=1900+dte.year;
        Gn._month=1+dte.mon;
        Gn._day=dte.mday;

        G(nping)=Gn;

        %Q
        Qn=struct;
          Qn.channel=nc;
          Qn.frequency=100e+3; %¿?
          Qn.transmitPower=1000; %¿?
          td=typecast(uint8(hd(41:44)'), "single");
          Qn.sampleInterval=0.3048*td/(soundSpeed*l);
          Qn.pulseLength=4*Qn.sampleInterval;
          Qn.soundVelocity=soundSpeed;
          Qn.absorptionCoefficient=0.0;
          Qn.temperature=temperature;
          Qn.count=n;
          dp=typecast(uint8(hd(61:64)'), "single"); %bad bottom detection in shallow areas
          Qn._depth=0.5*0.3048*dp;

        Q(nping)=Qn;

        ec=fread(f,l,"uint8");
        if( isempty(P) )
          P=nan(npings(1), l);
        endif
        P(nping,1:l)=(ec'-255)/log10(2);

      else
        fseek(f, l, "cof");
      endif

      ll=fread(f,1,'uint32');
      if( ll+4 ~= ftell(f) )
        error('Error in file')
      endif

    endfor
	%
  endfor

  P={P};
  Q={Q};

endfunction
