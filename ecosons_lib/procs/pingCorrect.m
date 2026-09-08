%
%PP=pingCorrect(P,d,d0,dmx, alpha,algn_dB)
%P: ping (dB)
%d: ping depth
%d0: reference depth
%dmx: depth the ping is stretched to
%alpha: attenuation coefficient (dB/meter); defaults to 0 (no attenuation)
%algn_dB: dB below the maximum to align the corrected ping to; defaults to 30 dB
%PP: corrected ping (dB)
function PP=pingCorrect(P,d,d0,dmx, alpha, algn_dB, l_bins)
 if( d>=d0 )

  PP=P+30*log10(d/d0);
  %!some improvement?
  %rr=(1/dmx)*[length(P):2*length(P)-1];
  %PP=P+20*log10(d/d0)+...
  %     10*log10((4+d*rr)./(4+d0*rr));

  if( exist('alpha') && length(alpha)==1 )
   PP=PP+2*alpha*d*(1+[0:length(P)-1]/length(P));
  endif

  if( ~exist('l_bins') )
    l_bins=4;
  endif

  if(dmx>d && d>d0)
   lp=l_bins*dmx/d;
   lp0=l_bins*dmx/d0;
   dlp=1+(lp0-lp);
   ckern=[ones(1,floor(dlp)), dlp-floor(dlp)];
   PP=10*log10( conv( power(10, 0.1*PP), ckern)/dlp ); %Use this one
   %!PP=20*log10( conv( power(10, 0.5*0.1*PP), ckern)/dlp ); %Test this one based on pressures?

   %ping alignment with the same 30 dB criterium
   if( ~exist('algn_dB') || length(algn_dB)==0 )
    algn_dB=30;
   endif

   [PPmx,mx_idx]=max(PP);
   p0=mx_idx;
   while(p0>1 && PPmx-PP(p0-1)<algn_dB)
    p0=p0-1;
   endwhile
    PP=[PP(p0:min(length(PP),length(P)+p0-1)),...
         nan(1,max(0,length(P)+p0-1-length(PP)))]; %should not overflow (but it might)
  endif

 else
  PP=nan(size(P));
 endif
endfunction
