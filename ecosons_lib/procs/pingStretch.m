%
%[PP1,PP2]=pingStretch(P,R,Rmx)
%P: ping (dB)
%R: ping depth (bins)
%Rmx: depth the ping is stretched to (bins)
%R2: apparent depth of the second echo (bins)
%PP1: stretched ping (dB)
%PP2: stretched second echo (dB)
function [PP,PP2]=pingStretch(P,R,Rmx, R2, l_bins)

%first and second bounces
PP=nan(1,Rmx);
PP2=nan(1,Rmx);

if( ~exist('l_bins') )
  l_bins=4;
  l_bins2=round(l_bins/2);
endif

if( R>4 )

 ll=min(R-l_bins2, size(P,2)-R);
 if(ll>1)

  PP=interp1([0:ll], P(R:R+ll), (R-l_bins2)/(Rmx-1)*[0:Rmx-1], 'linear', 'EXTRAP');
  %PP=interp1([0:ll], P(R:R+ll), (R-l_bins2)/(Rmx-1)*[0:Rmx-1], 'nearest', 'EXTRAP');

 endif

 if( exist('R2') && R2>R+4 )

  ll=min(R-2, size(P,2)-R2);
  if(ll>1)

   PP2=interp1([0:ll], P(R2:R2+ll), (R-2)/(Rmx-1)*[0:Rmx-1], 'linear', 'EXTRAP');

  endif

 endif

endif

endfunction
