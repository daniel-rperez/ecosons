%alpha=alphaAinslieMcColm(f,T,S,D,pH)
%Returns the sea water attenuation coefficient in dB/km
%Ref: Ainslie & McColm, J. Acoust. Soc. Am., Vol. 103, No. 3, March 1998
% f: frequency in Hz
% T: temperature in Celsius
% S: salinity in PSU
% D: depth (as a proxy to pressure)
% pH
% alpha: total absorption coefficient in dB/km
%%Version: 2026/09/10
function alpha=alphaAinslieMcColm(f,T,S,D,pH)
 f=f/1000; %f in kHz
 D=D/1000; %D in km
 %T in ªC
 %S in PSU ~= 35
 %pH ~= 8
 
 %Boric acid
 A1=0.106*exp((pH - 8)/0.56);
 P1=1;
 f1=0.78*sqrt(S/35).*exp(T/26);
 Boric=(A1.*P1.*f1.*power(f,2))./(power(f,2)+power(f1,2));

 %MgSO4 contribution
 A2=0.52*(S/35).*(1+T/43);
 P2=exp(-D/6);
 f2=42*exp(T/17);
 MgSO4=(A2.*P2.*f2*power(f,2))./(power(f,2)+power(f2,2));

 %Pure water contribution
 A3=0.00049*exp(-(T/27+D/17));
 P3=1;
 H2O=A3.*P3.*power(f,2);

 %Total absorption (dB/km)
 alpha=(Boric+MgSO4+H2O);
 alpha=alpha/1000; %dB/m

endfunction

