function [P,H,G]=fmt_lowranceSLG(fname)

  f=fopen(fname, 'rb', 'ieee-le');

  %read header
  slF=fread(f,1, 'uint16'); %format: 1 (SLG), 2 (SL2)
  slV=fread(f,1, 'uint16'); %version: 1, 2, 3
  %slB=fread(f,1, 'uint16');
  %sl_=fread(f,1, 'uint16');

  fclose(f);

  switch(slF)
    case 1
      %sl__=fread(f,1, 'int16');
      [P,H,G]=fmt_lowranceSL1(fname);
    case 2
      [P,H,G]=fmt_lowranceSL2(fname);
    otherwise
      error([ 'Format SL' num2str(slF) ' not yet supported'])
  endswitch

endfunction
