% function conf=loadConfig(fname)
% conf: 
% fname: 
%%Version: 2024/05/10
function conf=loadConfig(fname)

  %open file
  f=fopen(fname, "r");
  conf=struct();
  sub_k=[];
  nl=0;

  while( ~feof(f) )
    nl=nl+1;
    l=fgetl(f);
    l=strtrim(l);
    if( isempty(l) || l(1)==';' || l(1)=='#' || l(1)=='%' )
      continue
    endif

    %[s, ~, ~, ~, toks] = regexp (l, "^\\[([a-zA-Z0-9_]+)\\]$|^([a-zA-Z0-9_]+) *= *(true|false|on|off|t|f|'[^']+'|\"[^\"]+\"|[+-]*[0-9\\.]+( *, *[+-]*[0-9\\.]+)*|) *$");
    [s, ~, ~, ~, toks] = regexp (l, "^\\[('[^']+'|\"[^\"]+\"|[a-zA-Z][a-zA-Z0-9_]*)\\]$|^('[^']+'|\"[^\"]+\"|[a-zA-Z0-9_]+) *= *(true|false|on|off|t|f|@[a-zA-Z][a-zA-Z0-9_]*|'[^']+'|\"[^\"]+\"|[+-]*[0-9\\.]+[ei0-9+-]*( *, *[+-]*[0-9\\.]+[ei0-9+-]*)*|) *$", 'ignorecase');
    if( isempty(s) )
      error(["Unexpected syntax in " fname ':' num2str(nl)]);
    endif

    %sub-config?
    k=toks{1}{1};
    if( k(1)=="'" || k(1)=='"' ) %Extension: keys can be quoted strings (not limited to letters and numbers)
      k=k(2:end-1);
    endif
    if( length(toks{1}) == 1 )
      sub_k=k;
      continue
    endif

    %value string
    v=strtrim(toks{1}{2});

    %parse value string
    if( isempty(v) )
      vv=[];
    elseif( strcmp(tolower(v), 'true') || strcmp(tolower(v), 't') || strcmp(tolower(v), 'on') )
      vv=true;
    elseif( strcmp(tolower(v), 'false') || strcmp(tolower(v), 'f') || strcmp(tolower(v), 'off') )
      vv=false;
    elseif( (v(1)=='"' && v(end)=='"') || (v(1)=="'" && v(end)=="'") )
      vv=v(2:end-1);
    elseif( v(1)=='@' )
      vv=str2func(v);
    else
      vv=str2num(v);
    endif

    %append key-value pair to conf
    if( isempty(sub_k) )
      conf=setfield(conf, k, vv);
    else
      if( isfield(conf, sub_k) )
        cc=getfield(conf, sub_k);
      else
        cc=struct();
      endif
      cc=setfield(cc, k, vv);
      conf=setfield(conf, sub_k, cc);
    endif
  endwhile

  fclose(f);

endfunction

