% function saveConfig(conf, fname)
% conf: 
% fname: 
%%Version: 2024/05/10
function saveConfig(conf, fname)

  %open file
  if( ischar(fname) )
    f=fopen(fname, "w");
  else
    f=fname;
  endif

  for s=[false, true]
    for n=1:length(fieldnames(conf))

      k=fieldnames(conf){n};
        if( ~any(regexp(k, '^[a-zA-Z_][a-zA-Z_]+$')) )
          k=["'" k "'"];
        endif
      v=getfield(conf, k);
      
      if( s )
        if( isstruct(v) )
          if( ischar(fname) || fname == 1 )
            fprintf(f, '\n[%s]\n', k);
            saveConfig(v, f);
          else
            error(["Nested struct " k " not supported by config syntax!"]);
          endif
        endif
      else
        switch( class(v) )
        
          case 'char'
            if( index(v, '"') == 0 )
              fprintf(f, '%s="%s"\n', k, v);
            else
              fprintf(f, "%s='%s'\n", k, v);
            endif
          case 'double'
            if( length(v) == 1 )
              fprintf(f, '%s=%f\n', k, v);
            else
              fprintf(f, '%s=[%s]\n', k, num2str(v));
            endif
          case 'logical'
            if( length(v) == 1 )
              if( v )
                fprintf(f, '%s=true\n', k);
              else
                fprintf(f, '%s=false\n', k);
              endif
            else
              fprintf(f, '%s=[', k);
              for nv=1:length(nv)
                if( v(nv) )
                  fprintf(f, ' true');
                else
                  fprintf(f, ' false');
                endif
              endfor
              fprintf(f, ']\n');
            endif
          case 'function_handle'
            fprintf(f, '%s=@%s\n', k, func2str(v));
          case 'struct'
            fprintf(f, ';[%s] below\n', k);
          otherwise
            error(["Unrecognized type '" class(v) "' of " k " in config."]);
        endswitch
      endif
      
    endfor
  endfor

  if( ischar(fname) )
    fclose(f);
  endif

endfunction

