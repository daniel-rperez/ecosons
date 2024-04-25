%Reads a single-beam Simrad RAW sonar file into memory
% [P HS PS At Al]=fmt_simradRAW(fname, channel, mxp)
% P{}: channel cells with ping matrices: rows: ping no., columns: echo sample
% HS: channel cells with transducer headers
% PS: GPS + time data (time=-1: no data)
% At, Al: athwart- and along-ship electrical angles
% fname: RAW file filename
% RAW0: https://www.ngdc.noaa.gov/mgg/wcd/simradEK60manual.pdf
% RAW3: https://www.kongsberg.com/globalassets/maritime/km-products/product-documents/413763_ea640_ref.pdf
function [P,HS,PS,At,Al,W,filt]=fmt_simradRAW(fname)

 %Initialize returned values
 P={};
 W={};
 At={};
 Al={};
 HS={};
  lHS=struct;
  lChIds={};
 PS=struct;
  lPS=gpsRead("");

 %number of arguments
 request_As=[];
 request_Wf=[];

 %open file
 sonarFile=fopen(fname, 'rb');

 %sizes
 max_channels=0;
 max_counts=[0];

 %transceiver Z
 chTrcvZ=[0];
 filt=struct;

 pass=0;
 while(pass < 2)
  pass=pass+1;
  fseek(sonarFile, 0, "bof");

  %Initialize counter
  nPings=zeros(1,max_channels);

  while( ~feof(sonarFile) )

   dgrm=readRawDatagram(sonarFile);

   if( pass==1 )

    if( strcmp(dgrm.type, 'XML0') && isfield(dgrm.xml, 'Configuration') )
     tcvs=dgrm.xml.Configuration.Transceivers.Transceiver;
     for nt=1:length(tcvs)
      tchn=dgrm.xml.Configuration.Transceivers.Transceiver(nt).Channels.Channel;
      for nc=1:length(tchn)
       ch=find( strcmp(dgrm.xml.Configuration.Transceivers.Transceiver(nt).Channels.Channel(nc)._ChannelID, lChIds) );
       if( ~any(ch) )
        lChIds{end+1}=dgrm.xml.Configuration.Transceivers.Transceiver(nt).Channels.Channel(nc)._ChannelID;
        chTrcvZ(length(lChIds))=str2num(dgrm.xml.Configuration.Transceivers.Transceiver(nt)._Impedance);
       endif
      endfor
     endfor
    endif

%    if( strcmp(dgrm.type, 'XML0') && isfield(dgrm.xml, 'Channel') )
%    %if( strcmp(dgrm.type, 'XML0') && strcmp(dgrm.xml.__name, 'Channels') )
%     dgrm.xml
%     for nc=1:length(dgrm.xml.Channel)
%      ch=find( strcmp(dgrm.xml.Channel(nc)._ChannelID, lChIds) );
%      if( ~any(ch) )
%       lChIds{end+1}=dgrm.xml.Channel(nc)._ChannelID;
%      endif
%     endfor
%    endif

    if( strcmp(dgrm.type, 'RAW0') || strcmp(dgrm.type, 'RAW3') )

     %sample header
     lH=dgrm.sample;
     if( isfield(lH, 'channelId') )
       lH.channel=find( strcmp(lH.channelId, lChIds) );
     endif

     ch=lH.channel;

     if( ch > max_channels )
      max_counts(max_channels+1:ch)=0;
      max_channels=ch;
      nPings(ch)=0;
      request_Wf(ch)=0;
      request_As(ch)=false;
     endif

     if( isfield(dgrm, 'data') && isfield(dgrm.data, 'waveform') )
      if( nargout > 5 )
        request_Wf(ch)=size(dgrm.data.waveform,2);
      endif
      request_As(ch)=(nargout>3);
     endif

     if( isfield(dgrm, 'data') && isfield(dgrm.data, 'angleAthwartship') )
      request_As(ch)=(nargout>3);
     endif

     nPings(ch)=nPings(ch)+1;
     max_counts(ch)=max(max_counts(ch),lH.count);

    endif

   else

    switch(dgrm.type)

    case 'NME0'
     llPS=gpsRead(dgrm.nmea);
     if( llPS.time >=0 )
      lPS=llPS;
     endif

    case 'XML0'

     if( isfield(dgrm.xml, 'InitialParameter') )
      chs=dgrm.xml.InitialParameter.Channels;
      for nc=1:length(chs.Channel)
       ch=find( strcmp(chs.Channel(nc)._ChannelID, lChIds) );
       lHS(ch).channel=ch;
       lHS(ch).mode=str2num(chs.Channel(nc)._ChannelMode);
       lHS(ch).pulseForm=str2num(chs.Channel(nc)._PulseForm);
       if( isfield(chs.Channel(nc), '_Slope') )
        lHS(ch).slope=str2num(chs.Channel(nc)._Slope);
       else
        lHS(ch)=0.0;
       endif
       if( isfield(chs.Channel(nc), '_Frequency') && !isempty(chs.Channel(nc)._Frequency) )
        lHS(ch).frequency=str2num(chs.Channel(nc)._Frequency);
       elseif( isfield(chs.Channel(nc), '_FrequencyStart') )
        lHS(ch).frequency=[str2num(chs.Channel(nc)._FrequencyStart) str2num(chs.Channel(nc)._FrequencyEnd)];
       else
        lHS(ch).frequency=0;
        warning('Channel datagram without frequency?');
       endif
       lHS(ch).bandWidth=diff([lHS(ch).frequency 0])(1);
       lHS(ch).pulseLength=str2num(chs.Channel(nc)._PulseDuration);
       lHS(ch).sampleInterval=str2num(chs.Channel(nc)._SampleInterval);
       lHS(ch).transmitPower=str2num(chs.Channel(nc)._TransmitPower);
      endfor
     endif

     if( isfield(dgrm.xml, 'Parameter') )
      chn=dgrm.xml.Parameter.Channel;
      ch=find( strcmp(chn._ChannelID, lChIds) );
      lHS(ch).channel=ch;
      lHS(ch).mode=str2num(chn._ChannelMode);
      lHS(ch).pulseForm=str2num(chn._PulseForm);
      if( isfield(chn, '_Slope') )
       lHS(ch).slope=str2num(chn._Slope);
      else
       lHS(ch)=0.0;
      endif
      if( isfield(chn, '_Frequency') && !isempty(chn._Frequency) )
       lHS(ch).frequency=str2num(chn._Frequency);
      elseif( isfield(chn, '_FrequencyStart') )
       lHS(ch).frequency=[str2num(chn._FrequencyStart) str2num(chn._FrequencyEnd)];
      else
       lHS(ch).frequency=0;
       warning('Channel datagram without frequency?');
      endif
      lHS(ch).bandWidth=diff([lHS(ch).frequency 0])(1);
      lHS(ch).pulseLength=str2num(chn._PulseDuration);
      lHS(ch).sampleInterval=str2num(chn._SampleInterval);
      lHS(ch).transmitPower=str2num(chn._TransmitPower);
     endif

     if( isfield(dgrm.xml, 'Environment') ) %only one per file?
      for ch=1:length(lChIds)
       lHS(ch).transducerDepth=( depth=str2num(dgrm.xml.Environment._Depth) );
       lHS(ch).soundVelocity=str2num(dgrm.xml.Environment._SoundSpeed);
       lHS(ch).temperature=( temperature=str2num(dgrm.xml.Environment._Temperature) );
       lHS(ch).salinity=( salinity=str2num(dgrm.xml.Environment._Salinity) );
       lHS(ch).acidity=( acidity=str2num(dgrm.xml.Environment._Acidity) );
       %!!!disp(['f=' num2str(lHS(ch).frequency) ', T=' num2str(temperature) ', S=' num2str(salinity), ', D=' num2str(depth) ', pH=' num2str(acidity)])
       if( isfield(lHS(ch), 'frequency') )
        lHS(ch).absorptionCoefficient=alphaAinslieMcColm(lHS(ch).frequency,temperature,salinity,depth,acidity); %calculado a partir de la salinidad, el pH, la temperatura, ...
       else
        lHS(ch).absorptionCoefficient=0.0;
       endif
      endfor
     endif

    case 'MRU0'

     for ch=1:length(lChIds)
      lHS(ch).heave=dgrm.heave;
      lHS(ch).roll=dgrm.roll;
      lHS(ch).pitch=dgrm.pitch;
      lHS(ch).heading=dgrm.heading;
     endfor

    case 'CON0'

     if( length(dgrm.transducer) > 0 )
      for ch=1:length(lChIds)
       lHS(ch).beamType=dgrm.transducer(1).beamType;
       lHS(ch).gain=dgrm.transducer(1).gain;
       lHS(ch).equivalentBeamAngle=dgrm.transducer(1).equivalentBeamAngle;
      endfor
     endif

    case 'RAW0'

     %sample header
     lH=dgrm.sample;
     ch=lH.channel;

     nPings(ch)=nPings(ch)+1;

     %header and last GPS data
     HS{ch}(nPings(ch))=lH;
     PS(nPings(ch))=lPS;

     %sample power
     if( isfield(dgrm.data, "power") )

      ll=min(lH.count, length(dgrm.data.power));
      P{ch}(nPings(ch),1:ll)=dgrm.data.power(1:ll);

     endif

     %sample angles
     if( request_As(ch) ...
      && isfield(dgrm.data, "angleAthwartship") ...
      && isfield(dgrm.data, "angleAlongship") )

      ll=min(lH.count, length(dgrm.data.power));
      At{ch}(nPings(ch),1:ll)=dgrm.data.angleAthwartship(1:ll);
      Al{ch}(nPings(ch),1:ll)=dgrm.data.angleAlongship(1:ll);

     endif

    case 'RAW3'

     %sample header
     lH=dgrm.sample;
     lH.channel=find( strcmp(lH.channelId, lChIds) );
     ch=lH.channel;
     ff=fieldnames(lHS);
     %!!!disp([ "'" lH.channelId "'" " ?= " "'" lChIds{1} "'" " | " "'" lChIds{2} "'" " : " num2str(ch) ])
     for nf=1:length(ff)
      lH=setfield(lH, ff{nf}, getfield(lHS(ch),ff{nf}));
     endfor

     nPings(ch)=nPings(ch)+1;

     %header and last GPS data
     HS{ch}(nPings(ch))=lH;
     PS(nPings(ch))=lPS;

     %sample power
     if( isfield(dgrm.data, "power") )

      ll=min(lH.count, length(dgrm.data.power));
      P{ch}(nPings(ch),1:ll)=dgrm.data.power(1:ll);

      %sample angles
      if( request_As(ch) ...
       && isfield(dgrm.data, "angleAthwartship") ...
       && isfield(dgrm.data, "angleAlongship") )

       ll=min(lH.count, length(dgrm.data.power));
       At{ch}(nPings(ch),1:ll)=dgrm.data.angleAthwartship(1:ll);
       Al{ch}(nPings(ch),1:ll)=dgrm.data.angleAlongship(1:ll);

      endif

     endif

     %waveform
     if( isfield(dgrm.data, "waveform") )

      ll=min(lH.count, size(dgrm.data.waveform,1));
      W{ch}(nPings(ch),1:ll,:)=(1000 + 75)/1000 * dgrm.data.waveform(1:ll,:);

      dTX=round(2*lH.pulseLength/lH.sampleInterval); %!número de bins que contienen la señal transmitida
      P{ch}(nPings(ch),1:ll-dTX)=10*log10( 0.5*power(abs( mean(dgrm.data.waveform(1+dTX:ll,:),2) ), 2)/75 );

      if( request_As(ch) )
       if( isfield(lH, 'beamType') )
        bt=lH.beamType;
       else
        bt=0*(request_Wf(ch)==1)+1*(request_Wf(ch)==4)+17*(request_Wf(ch)==3);
       endif
       if( bt>0 )
        [paAt,paAl]=phaseAngle(squeeze(W{ch}(nPings(ch),1:ll,:)), bt);
        At{ch}(nPings(ch),1:ll-dTX)=paAt(1+dTX:end);
        Al{ch}(nPings(ch),1:ll-dTX)=paAl(1+dTX:end);
       else
        At{ch}(nPings(ch),1:ll-dTX)=nan(1,ll-dTX);
        Al{ch}(nPings(ch),1:ll-dTX)=nan(1,ll-dTX);
       endif
      endif
     endif

    case 'FIL1'
     ch=find( strcmp(dgrm.channelID, lChIds) );
     filt(ch).noOfCoefficients(dgrm.stage)=dgrm.noOfCoefficients;
     filt(ch).decimationFactor(dgrm.stage)=dgrm.decimationFactor;

    endswitch

   endif

  endwhile

  if( pass == 1 )

   for ch=1:max_channels
    if( request_Wf(ch) )
     P{ch} = nan(nPings(ch),max_counts(ch));
     W{ch} = nan(nPings(ch),max_counts(ch), request_Wf(ch));
    else
     P{ch} = nan(nPings(ch),max_counts(ch));
    endif
    if( request_As(ch) )
     Al{ch}=nan(nPings(ch),max_counts(ch));
     At{ch}=nan(nPings(ch),max_counts(ch));
    endif
   endfor

  endif

 endwhile %pass
 fclose(sonarFile);

endfunction


function lPS=gpsRead(nmea)
 lPS.time=-1.0;
 lPS.latitude=0.0;
 lPS.longitude=0.0;

 if( length(nmea) > 5 )
  ss=strsplit(nmea(4:end), ',', 'collapsedelimiters', false);
  switch(ss{1})

  %https://gpsd.gitlab.io/gpsd/NMEA.html#_gga_global_positioning_system_fix_data

  case 'GGA' %GPS fix data: UTC lat lon
   %check for format errors (missing commas, non-numeric values) or missing coordinates
   if( (length(ss) == 15) && ...
       ~isempty(ss{2}) && ~isempty(ss{3}) && ~isempty(ss{5}) && ...
       ~isalpha(ss{2}(1)) && ~isalpha(ss{3}(1)) && ~isalpha(ss{5}(1)) )

    ll=ss{2};
    lPS.time=str2num(ll(1:2))+str2num(ll(3:4))/60+str2num(ll(5:6))/3600;

    ll=ss{3};
    lPS.latitude=str2num(ll(1:2))+str2num(ll(3:end))/60;
     if( ss{4} == 'S' )
      lPS.latitude=-lPS.latitude;
     endif
    ll=ss{5};
    lPS.longitude=str2num(ll(1:3))+str2num(ll(4:end))/60;
     if( ss{6} == 'W' )
      lPS.longitude=-lPS.longitude;
     endif
   endif

  case 'GLL' %Geographic position: latitude longitude UTC
   ll=ss{2};
   lPS.latitude=str2num(ll(1:2))+str2num(ll(3:end))/60;
    if( ss{3} == 'S' )
     lPS.latitude=-lPS.latitude;
    endif
   ll=ss{4};
   lPS.longitude=str2num(ll(1:3))+str2num(ll(4:end))/60;
    if( ss{5}(1) == 'W' )
     lPS.longitude=-lPS.longitude;
    endif
   if(length(ss)>6)
     ll=ss{6};
     lPS.time=str2num(ll(1:2))+str2num(ll(3:4))/60+str2num(ll(5:6))/3600;
   else
     lPS.time=-1;
   endif

  case 'GNS' %Fix data: UTC latitude longitude
   ll=ss{2};
   lPS.time=str2num(ll(1:2))+str2num(ll(3:4))/60+str2num(ll(5:6))/3600;

   ll=ss{3};
   lPS.latitude=str2num(ll(1:2))+str2num(ll(3:end))/60;
    if( ss{4} == 'S' )
     lPS.latitude=-lPS.latitude;
    endif
   ll=ss{5};
   lPS.longitude=str2num(ll(1:3))+str2num(ll(4:end))/60;
    if( ss{6} == 'W' )
     lPS.longitude=-lPS.longitude;
    endif

  case 'GXA' %Transit position: UTC latitude longitude
   ll=ss{2};
   lPS.time=str2num(ll(1:2))+str2num(ll(3:4))/60+str2num(ll(5:6))/3600;

   ll=ss{3};
   lPS.latitude=str2num(ll(1:2))+str2num(ll(3:end))/60;
    if( ss{4} == 'S' )
     lPS.latitude=-lPS.latitude;
    endif
   ll=ss{5};
   lPS.longitude=str2num(ll(1:3))+str2num(ll(4:end))/60;
    if( ss{6} == 'W' )
     lPS.longitude=-lPS.longitude;
    endif

  endswitch

 endif

endfunction


function dgrm=readRawDatagram(sonarFile, ctx)
 dgrm=struct;
 dgrm.length=fread(sonarFile,1,'int32','ieee-le');
  if( length(dgrm.length) == 0 || dgrm.length==0 )
   dgrm.type='END0';
   return;
  endif
 dgrm.type=sprintf('%c', fread(sonarFile,4,'uchar'));
 dgrm.date(2)=fread(sonarFile,1,'int32','ieee-le');
 dgrm.date(1)=fread(sonarFile,1,'int32','ieee-le');

 switch(dgrm.type)
 case 'CON0'
  hdr=struct;
  hdr.surveyName=sprintf('%c', fread(sonarFile,128,'uchar'));
  hdr.transectName=sprintf('%c', fread(sonarFile,128,'uchar'));
  hdr.sounderName=sprintf('%c', fread(sonarFile,128,'uchar'));
  hdr.version=sprintf('%c', fread(sonarFile, 30,'uchar'));
   fread(sonarFile, 98,'uchar');
  hdr.transducerCount=fread(sonarFile, 1,'int32','ieee-le');

  dgrm.header=hdr;

  dgrm.transducer=struct;
  for ntransducer=1:dgrm.header.transducerCount
   tr_des=struct;
   tr_des.channelId=sprintf('%c', fread(sonarFile,128,'uchar'));
   tr_des.beamType=fread(sonarFile, 1,'int32','ieee-le');
   tr_des.frequency=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.gain=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.equivalentBeamAngle=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.beamAlongship=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.beamAthwartship=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.sensitivityAlongship=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.sensitivityAthwartship=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.offsetAlongship=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.offsetAthwartship=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.posX=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.posY=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.posZ=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.dirX=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.dirY=fread(sonarFile, 1,'float32','ieee-le');
   tr_des.dirZ=fread(sonarFile, 1,'float32','ieee-le');
    fread(sonarFile, 8,'uchar');
   tr_des.pulseLengthTable=fread(sonarFile, 5,'float32','ieee-le');
    fread(sonarFile, 8,'uchar');
   tr_des.gainTable=fread(sonarFile, 5,'float32','ieee-le');
    fread(sonarFile, 8,'uchar');
   tr_des.saCorrectionTable=fread(sonarFile, 5,'float32','ieee-le');
    fread(sonarFile, 8,'uchar');
   tr_des.gptSofwareVersion=sprintf('%c', fread(sonarFile,16,'uchar'));
    fread(sonarFile, 20,'uchar');


   dgrm.transducer(ntransducer)=tr_des;
  endfor

 case 'NME0'
  dgrm.nmea=sprintf('%c',fread(sonarFile,dgrm.length-3*4,'uchar'));
  %!!!disp(['NMEA=' dgrm.nmea]);

 case 'XML0'
  txt=sprintf('%c',fread(sonarFile,dgrm.length-3*4,'uchar'));
  try
    dgrm.xml=parseXML(txt);
  catch
    disp(['XML: ' txt]);
  end_try_catch
  %XML0: Transducers.Transducer(1:N), Transceivers.Transceiver(1:N), etc.
  %XML0: Channel(1:N).{Frequency,PulseDuration,SampleInterval,TransmitPower,SoundVelocity} también un único Channel
  %XML0: Environment.{Depth, SoundSpeed}
  %XML0: Sensor
  %XML0: Ping

 case 'TAG0'
  dgrm.text=sprintf('%c',fread(sonarFile,dgrm.length-3*4,'uchar'));
   if( length(dgrm.text)>80 )
    dgrm.text=dgrm.text(1:80);
   endif
   %!!!disp(['TAG0=' dgrm.text]);

 case 'MRU0'
  dgrm.heave=fread(sonarFile,1,'float32');
  dgrm.roll=fread(sonarFile,1,'float32');
  dgrm.pitch=fread(sonarFile,1,'float32');
  dgrm.heading=fread(sonarFile,1,'float32');

 case 'RAW0'
  samp=struct;

  samp.channel=fread(sonarFile, 1,'int16','ieee-le');
  samp.mode=fread(sonarFile, 1,'int16','ieee-le');
  samp.transducerDepth=fread(sonarFile, 1,'float32','ieee-le');
  samp.frequency=fread(sonarFile, 1,'float32','ieee-le');
  samp.transmitPower=fread(sonarFile, 1,'float32','ieee-le');
  samp.pulseLength=fread(sonarFile, 1,'float32','ieee-le');
  samp.bandWidth=fread(sonarFile, 1,'float32','ieee-le');
  samp.sampleInterval=fread(sonarFile, 1,'float32','ieee-le');
  samp.soundVelocity=fread(sonarFile, 1,'float32','ieee-le');
  samp.absorptionCoefficient=fread(sonarFile, 1,'float32','ieee-le');
  samp.heave=fread(sonarFile, 1,'float32','ieee-le');
  samp.roll=fread(sonarFile, 1,'float32','ieee-le');
  samp.pitch=fread(sonarFile, 1,'float32','ieee-le');
  samp.temperature=fread(sonarFile, 1,'float32','ieee-le');
   fread(sonarFile, 12,'uchar');
  samp.offset=fread(sonarFile, 1,'int32','ieee-le');
  samp.count=fread(sonarFile, 1,'int32','ieee-le');

  dat=struct;
  if( samp.mode==1 || samp.mode==3 )
   x=fread(sonarFile, samp.count,'int16','ieee-le');
   dat.power=10*log10(2)*x/256;
  endif
  if( samp.mode==2 || samp.mode==3 )
   x=fread(sonarFile, 2*samp.count,'int8','ieee-le');
   dat.angleAthwartship=180*x(1:2:end)/128; %most significant byte
   dat.angleAlongship=180*x(2:2:end)/128;   %least significant byte
  endif

  dgrm.sample=samp;
  dgrm.data=dat;

 case 'RAW3'
  samp.channelId=sprintf('%c', fread(sonarFile,128,'uchar'));
   samp.channelId=samp.channelId( samp.channelId ~= 0 );
  samp.dataType=fread(sonarFile, 1,'int16','ieee-le');
  samp.spare=fread(sonarFile, 1,'int16','ieee-le');
  samp.offset=fread(sonarFile, 1,'int32','ieee-le');
  samp.count=fread(sonarFile, 1,'int32','ieee-le');

  dat=struct;
  if( bitand(samp.dataType,1) ) %Power: integer "compressed"
   x=fread(sonarFile, samp.count,'int16','ieee-le');
   dat.power=10*log10(2)*x/256;
  endif
  if( bitand(samp.dataType,2) ) %Angle: integer encoded
   x=fread(sonarFile, 2*samp.count,'int8','ieee-le');
   dat.angleAthwartship=180*x(1:2:end)/128; %most significant byte
   dat.angleAlongship=180*x(2:2:end)/128;   %least significant byte
  endif
  if( bitand(samp.dataType,4) ) %Waveform: complex float 16
   cperSample=floor(samp.dataType/256);
   x=fread(sonarFile, 2*cperSample*samp.count,'float16','ieee-le');
   dat.waveform=reshape(x(1:2:end-1)+j*x(2:2:end), [cperSample, samp.count])';
  endif
  if( bitand(samp.dataType,8) ) %Waveform: complex float 32
   cperSample=floor(samp.dataType/256);
   x=fread(sonarFile, 2*cperSample*samp.count,'float32','ieee-le');
   dat.waveform=reshape(x(1:2:end-1)+j*x(2:2:end), [cperSample, samp.count])';
  endif

  dgrm.sample=samp;
  dgrm.data=dat;

case 'FIL1'
  dgrm.stage=fread(sonarFile, 1,'int16','ieee-le');
  fread(sonarFile, 1,'char'); %!!!spec says 2 chars!
  dgrm.filterType=fread(sonarFile, 1,'char');
  dgrm.channelID=char( fread(sonarFile, 128,'char')' );
   dgrm.channelID=dgrm.channelID( dgrm.channelID ~= 0 );
  dgrm.noOfCoefficients=fread(sonarFile, 1,'int16','ieee-le');
  dgrm.decimationFactor=fread(sonarFile, 1,'int16','ieee-le');
  x=fread(sonarFile, 2*dgrm.noOfCoefficients,'float32','ieee-le');
  dgrm.coefficients=x(1:2:end-1)+j*x(2:2:end);

 otherwise
  disp(['Unknown datagram: ' dgrm.type]);
  fread(sonarFile,dgrm.length-3*4,'uchar');
 endswitch

 dgrm_lc=fread(sonarFile,1,'int32','ieee-le');
 if( dgrm_lc ~= dgrm.length )
  error('File does not conform Simrad Raw Format');
 endif

endfunction

%Sound attenuation
function alpha=alphaAinslieMcColm(f,T,S,D,pH) %Ainslie & McColm, J. Acoust. Soc. Am., Vol. 103, No. 3, March 1998
 f=f/1000; %f in kHz
 D=D/1000; %D in km
 %T in ªC
 %S in PSU ~= 35
 %pH ~= 8

 %Boric acid
 A1=0.106*exp((pH - 8)/0.56);
 P1=1;
 f1=0.78*sqrt(S/35)*exp(T/26);
 Boric=(A1*P1*f1*power(f,2))./(power(f,2)+f1^2);

 %MgSO4 contribution
 A2=0.52*(S/35)*(1+T/43);
 P2=exp(-D/6);
 f2=42*exp(T/17);
 MgSO4=(A2*P2*f2*power(f,2))./(power(f,2)+f2^2);

 %Pure water contribution
 A3=0.00049*exp(-(T/27+D/17));
 P3=1;
 H2O=A3*P3*power(f,2);

 %Total absorption (dB/km)
 alpha=(Boric+MgSO4+H2O);

endfunction

%phase angle from Waveform for different beam types
% currently only beamtype 1 with 4 signals
function [paAt,paAl]=phaseAngle(W, bt)
 sW=size(W);
 L=prod(sW(1:end-1));
 paAt=reshape(zeros(L,1), [1,sW(1:end-1)]);
 paAl=reshape(zeros(L,1), [1,sW(1:end-1)]);
 switch( bt )
  case 1 %Transducers having four sectors: Starboard Aft, Port Aft, Port Fore, Starboard Fore
    paAl(:)=arg( conj(W(:,1)+W(:,2)).*(W(:,3)+W(:,4)) );
    paAt(:)=arg( conj(W(:,2)+W(:,3)).*(W(:,1)+W(:,4)) );
  case 17 %Transducers having three sectors: Starboard Aft, Port Aft, Forward
    paAl(:)=( arg( conj(W(:,1)).*W(:,3) ) + arg( conj(W(:,2)).*W(:,3) ) )/sqrt(3);
    paAt(:)=( arg( conj(W(:,2)).*W(:,3) ) - arg( conj(W(:,1)).*W(:,3) ) );
  case {49,65,81} %Transducers having three sectors and a centre element: Starboard Aft, Port Aft, Forward, Centre
    paAl(:)=( arg( conj(W(:,1)+W(:,4)).*(W(:,3)+W(:,4)) ) + arg( conj(W(:,2)+W(:,4)).*(W(:,3)+W(:,4)) ) )/sqrt(3);
    paAt(:)=( arg( conj(W(:,2)+W(:,4)).*(W(:,3)+W(:,4)) ) - arg( conj(W(:,1)+W(:,4)).*(W(:,3)+W(:,4)) ) );
  case 97 %Transducers having four sectors: Fore Starboard, Aft Port, Aft starboard, Fore Port
    paAt(:)=arg( conj(W(:,2)).*W(:,1) );
    paAl(:)=arg( conj(W(:,4)).*W(:,3) );
 endswitch
endfunction
