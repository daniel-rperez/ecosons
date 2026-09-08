%Gets the first hits of an array of pings in a
% per ping bases (no ping comparison is performed)
%[hit]=getFirstHit(P,Q, ndB, nndB)
% P: matrix holding one ping per row
% Q: acquisition data corresponding to the pings in P
% nearF: near field in meters
% ndB: boundary limit of the main lobe (below 5 m depth)
% nndB: background noise limit
% Note: ndB~30dB, according to the directivity function
%  main/secondary lobe criterium; nndB~60dB (¿?)
%%Version: 2024/02/24
function [hit]=getFirstHit(P,Q,nearF, ndB, nndB)

 %near field index
 knf=2*floor( nearF/(Q(1).soundVelocity*Q(1).sampleInterval) );
 pP=P(:,knf:size(P,2));
 [maxP,kmax]=max(pP');
 hit=kmax;
 for p=1:size(kmax,2)
  for k=kmax(p):-1:1
   if( pP(p,k) >= maxP(p)-ndB )
    hit(p)=k;
   elseif( pP(p,k) < maxP(p)-nndB )
    break;
   endif
  endfor
 endfor

 hit=hit+knf-1; %-1 to correct knf

endfunction

