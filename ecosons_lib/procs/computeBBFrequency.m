%[fW, f_sel, P, Pmf]=computeBBFrequency(HS, W, f_rep, useWtemplate, withCorrections, calXML)
%
%
function [fW, f_sel, P, Pmf]=computeBBFrequency(HS, W, f_rep, useWtemplate, withCorrections, calXML)

  pkg load signal

  %Pulse length (in seconds and in bins)
  dts=HS(1).sampleInterval;
  dT=HS(1).pulseLength;
  ndT=round(dT / dts );
  %!ndT=2*round(1 + dT / (2*dts) ); %generous length

  if( ~exist('useWtemplate') )
    useWtemplate=false;
  else
    if( ~useWtemplate )
      error("Use of simulated template is not implemented")
    endif
  endif

  %Find initial chirp
  hchirp0=mean(W(1,1:3*ndT,:), 3);
  pchirp0=hchirp0.*conj(hchirp0);

  %Find center of chirp
  n_c=floor( sum( pchirp0 .* [1:length(pchirp0)] )/sum(pchirp0) );
  rchirp=floor(n_c-ndT):floor(n_c+ndT);
  rchirp=rchirp(0<rchirp);

  hchirp=mean(W(1,rchirp,:), 3);
  hchirp_acr=conv(hchirp, conj(hchirp(end:-1:1)), 'same');

  pow_scale=norm(hchirp)/norm(hchirp_acr); %to ensure energy conservation
  hchirp_acr=pow_scale*hchirp_acr;

  %Maximum of the chirp transform: set as t=0 (before first sample)
  [~,n_cc]=max(abs(conv(hchirp0, conj(hchirp(end:-1:1)), 'same')));

  %Time base
  t_chirp=[1:length(hchirp)]*HS(1).sampleInterval;
    t_chirp=t_chirp - mean(t_chirp);

  %Slipt beam sector average
  hsignal=mean(W, 3);

  %Matched filter
  hsignal_uconv=conv2(hsignal, conj(hchirp(end:-1:1)), 'same');
    hsignal_uconv=pow_scale*hsignal_uconv;

  %Hann window
  wHann=hann(length(hchirp))'; %as row vector

  %Frequency selection
  if( length(f_rep)==0 )
    f_sel=[ HS(1).frequency(1) : 1000 : HS(1).frequency(2) ];
  else
    f_sel=f_rep( (HS(1).frequency(1) <= f_rep)  &  (f_rep <= HS(1).frequency(2)) );
  endif

  if( exist('withCorrections') )
    if( isbool(withCorrections) )
      withCorrections=40;
    endif
  else
    withCorrections=0;
  endif

  %Calibration
  if( exist('calXML') )
    f_mean=mean(HS(1).frequency);
    cw=HS(1).soundVelocity;
    l_mean=cw/f_mean;
    pwEm=HS(1).transmitPower;
    Ffc=str2num( strrep(callXML.Root.Calibration.CalibrationResults.Frequency.__content, ';', ',') );
    Gfc=str2num( strrep(callXML.Root.Calibration.CalibrationResults.Gain.__content, ';', ',') );
    Zfc=str2num( strrep(callXML.Root.Calibration.CalibrationResults.Impedance.__content, ';', ',') );
    EqBA=10*log10(2 ...
           * str2num( strrep(callXML.Root.Calibration.CalibrationResults.BeamWidthAlongship.__content, ';', ',') ) ...
           * str2num( strrep(callXML.Root.Calibration.CalibrationResults.BeamWidthAthwartship.__content, ';', ',') ) ...
           );
    r_load=str2num( xmlCal.Root.Calibration.Common.Transceiver.Impedance.__content );

    if( min(Ffc)<=f_mean && f_mean<=max(Ffc) )
      g_sel=power(10, interp1(Ffc, Gfc, f_sel, 'spline', 'extrap')/10);
      z_sel=interp1(Ffc, Zfc, f_sel, 'spline', 'extrap');
      ea_sel=interp1(Ffc, EqBA, f_sel, 'spline', 'extrap');
      g_mean=power(10, mean(Gfc)/10);
      z_mean=mean(Zfc);
      ea_mean=mean(ea_sel);
    else
      g_sel=power(10, 20/10)*ones(size(f_sel));
      g_mean=power(10, 20/10);
      z_sel=75*ones(size(f_sel))
      z_mean=75;
      ea_sel=10*log10(2*17.3^2*(38000 ./ f_sel));
      ea_mean=mean(ea_sel);
    endif

  else
    f_mean=mean(HS(1).frequency);
    cw=HS(1).soundVelocity;
    l_mean=cw/f_mean;
    pwEm=HS(1).transmitPower;
    g_sel=power(10, 20/10)*ones(size(f_sel));
    g_mean=power(10, 20/10);
    r_load=5400; %!!!
    z_sel=75*ones(size(f_sel));
    z_mean=75;
    ea_sel=10*log10(2*17.3^2*(38000 ./ f_sel)); %!an approximation based on 38 kHz calibration
    ea_mean=mean(ea_sel);
  endif

  %Frequency cube
  fW=nan([length(f_sel), size(hsignal)]);
  fC=nan(size(f_sel));
  for nf=1:length(f_sel)
    w_f=wHann.*exp(-2i*pi*f_sel(nf)*t_chirp); %conjugated Fourier-Hann window!
    fC(nf)=hchirp_acr*w_f'; %Fourier component of the chirp autocorrelation (mind the traspose as it also conjugates)
    fWnf=conv2(hsignal_uconv,w_f, 'same')/fC(nf); %Normalized Fourier transform of the matched filtered signal
    fW(nf,:)=(r_load+z_sel(nf))/((1000+75)*r_load) * fWnf(:); %Correct electrical load gain and apply calibration (mechanical) gain
  endfor

  %Cut and reflect
  fW=fW(:,:,n_cc+1:end);

  %Power by frequency
  if( nargout>=3 )
    %Acoustic power received
    pW=nan(size(fW));
    for nf=1:length(f_sel)
      pW(nf,:)=0.5*size(W,3)*power(abs(fW(nf,:)),2)/( z_sel(nf) * g_sel(nf)^2 );
    endfor

    P={}; %convert to dB and separate by frequency
    for nf=1:length(f_sel)
      P{nf}=10*log10( reshape(pW(nf,:,:), [size(pW,2), size(pW,3)]) ) ...
           - 10*log10(pwEm*power(cw/(4*pi*f_sel(nf)),2)) ...
           + withCorrections*log10(0.5*[1:size(pW,3)]*cw*dts) ...
           - ea_sel(nf);
    endfor

    clear pW
  endif

  if( nargout>=4 )
    Pmf=10*log10( 0.5*power(abs(hsignal_uconv(:,n_cc+1:end)),2)/( z_mean * g_mean^2 ) ) ...
       - 10*log10(pwEm*power(cw/(4*pi*f_mean),2)) ...
       + withCorrections*log10(0.5*[1:size(fW,3)]*cw*dts) ...
       - ea_mean;
  endif

endfunction

