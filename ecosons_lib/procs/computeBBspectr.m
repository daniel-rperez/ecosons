function [fW, f_sel, P, Pmf, hchirp_acr]=computeBBspectr(HS, W, f_sel)

  W=sum(W,3);
  fW=W;

  dts=HS(1).sampleInterval;
  dT=HS(1).pulseLength;
  ndT=3*round(1 + dT / dts ); %generous length

  wwp=W(1:ndT);

  [wwp_mx,wwp_px]=max(abs(wwp));
  n1=wwp_px-1; while(n1>1 && abs(wwp(n1))>wwp_mx/256), n1=n1-1; endwhile
  n2=wwp_px+1; while(n2<=length(wwp) && abs(wwp(n2))>wwp_mx/256), n2=n2+1; endwhile
  wwp=wwp(n1:n2);
  tp=[0:length(wwp)-1]*HS.sampleInterval;

  if( ~exist('f_sel') )
    f_sel=[HS.frequency(1):1000:HS.frequency(2)];
  endif

  P={};
  for n=1:length(f_sel)
    ex=exp(1i*f_sel(n)*tp) .* hanning(length(tp))';
      ex=ex/norm(ex);

    fwwpn=sum(wwp.*ex);
    x=conv(W,ex(end:-1:1));
    [~,b]=max(abs(x));
    x=x(b:end);
    if(length(x)<length(W))
      x=[x zeros(1,length(W)-length(x))];
    endif
    x=x(1:length(W));

    P{n}=20*log10(abs(x/fwwpn));

  endfor

  hchirp_acr=conv(wwp,conj(wwp(end:-1:1)));
  Pmf=20*log10(abs(W));

endfunction

