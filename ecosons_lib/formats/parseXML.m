% [xmlout, tagC, tagE]=parseXML(txt, xmlin)
% Called as:
%  xmlout = parseXML(txt)
% Parses a text containing an xml document
% txt: xml document or fragment (complete elements)
% xmlin: the tag holding the document (or an empty document)
%  used only internally for recursive calls
% xmlout: octave struct holding the tag/document
%  xmlout.__name: tag name
%  xmlout.__prefix: tag prefix
%  xmlout.__attribute: attribute
%  xmlout.child: child xml tag (field name is tag name)
%  xmlout.children(): same tag children list (field names are tag names)
% tagC: boolean indicating whether the tag has been closed
%  used only internally for recursive calls
% tagE: position of txt where parsing has stopped (next character)
%  used only internally for recursive calls
function [xmlout, tagC, tagE]=parseXML(txt, xmlin, tagB)
  
  %straighten the text out (no newlines)
  if( ~exist('xmlin') )
    txt=strrep(strrep(txt, "\n", " "), "\r", "");

    xmlout=struct;
    xmlout.__name='.';
  else
    xmlout=xmlin;
  endif
    
  tagC=false;
  
  if( exist('tagB') )
    p=tagB;
  else
    p=1;
  endif

  while( p < length(txt) )
  
    %!
    %disp( [num2str(p) ' / ' num2str(length(txt))] )
    
    stxt=txt(p:min([p+1024, length(txt)]));
    try
      regexp(stxt, '^.*$'); %!fails if stxt does not include complete utf-8 characters
    catch
      stxt=stxt(1:end-1);
    end_try_catch
    %closing tag: return
    %[S,E]=regexp(txt(p:end), ['^</' xmlout.__name '\s*>'], 'once');
    [S,E]=regexp(stxt, ['^</' xmlout.__name '\s*>'], 'once');
    if( any(S) )
      tagC=true;
      p=p+E;
      break;
    endif

    %spaces: ignore
    %[S,E]=regexp(txt(p:end), '^\s+', 'once');
    [S,E]=regexp(stxt, '^\s+', 'once');
    if( any(S) )
      p=p+E;
      continue;
    endif
    
    %initial XML declaration: ignore
    %[S,E]=regexp(txt(p:end), '^<\?xml\s[^?]+\?>', 'once');
    [S,E]=regexp(stxt, '^<\?xml\s[^?]+\?>', 'once');
    if( any(S) )
      p=p+E;
      continue;
    endif
    
    %tag with attributes
    %[S,E]=regexp(txt(p:end), ['^<([a-zA-Z0-9:]+)\s*(([a-zA-Z0-9:]+)=("([^"]+)"|' "'" '([^' "'" ']+)' "'" ')\s*)+\s*/*>'], 'once');
    [S,E]=regexp(stxt, ['^<([a-zA-Z0-9:]+)\s*(([a-zA-Z0-9:]+)=("([^"]+)"|' "'" '([^' "'" ']+)' "'" ')\s*)+\s*/*>'], 'once');
    if( any(S) )

      %[sxml,stagC,stagE]=parseXMLTag(txt(p:p+E-1));
      [sxml,stagC,stagE]=parseXMLTag(stxt(1:E));
      
      p=p+E;
      
      %if the tag was not autoclosed read the contents
      if( ~stagC )
        [sxml, stagC, stagE]=parseXML(txt, sxml, p);
        p=stagE;
      endif
      
      %append field (if there was a previous one, add as a struct array)
      if( isfield(xmlout, sxml.__name) )
        f=getfield(xmlout, sxml.__name);
        lf=length(f);
        fn=fieldnames(sxml);
        for nf=1:length(fn)
          f=setfield(f, {lf+1}, fn{nf}, getfield(sxml,fn{nf}));
        endfor
        %f(end+1)=sxml;
        xmlout=setfield(xmlout, sxml.__name, f);
      else
        xmlout=setfield(xmlout, sxml.__name, sxml);
      endif
      
      continue
    endif

    %tag without attributes
    %[S,E]=regexp(txt(p:end), ['^<([a-zA-Z0-9:]+)\s*/*>'], 'once');
    [S,E]=regexp(stxt, ['^<([a-zA-Z0-9:]+)\s*/*>'], 'once');
    if( any(S) )

      %[sxml,stagC,stagE]=parseXMLTag(txt(p:p+E-1));
      [sxml,stagC,stagE]=parseXMLTag(stxt(1:E));
      
      p=p+stagE;
      
      %if the tag was not autoclosed read the contents
      if( ~stagC )
        [sxml, stagC, stagE]=parseXML(txt, sxml, p);
        p=stagE;
      endif
      
      %append field (if there was a previous one, add as a struct array)
      if( isfield(xmlout, sxml.__name) )
        f=getfield(xmlout, sxml.__name);
        lf=length(f);
        fn=fieldnames(sxml);
        for nf=1:length(fn)
          f=setfield(f, {lf+1}, fn{nf}, getfield(sxml,fn{nf}));
        endfor
        %f(end+1)=sxml;
        xmlout=setfield(xmlout, sxml.__name, f);
      else
        xmlout=setfield(xmlout, sxml.__name, sxml);
      endif
      
      continue
    endif
    
    %text in between: ignore?
    [S,E]=regexp(txt(p:end), '^[^<]+', 'once');
    if( any(S) )
      %! do not ignore
      if( isfield(xmlout, '__content') )
        xmlout.__content=[xmlout.__content ' ' strtrim(txt(p:p+E-1))];
      else
        xmlout.__content=strtrim(txt(p:p+E-1));
      endif
      p=p+E;
      continue
    endif
    
    %spare <: ignore error
    if( txt(p) == '<' )
     p=p+1;
    endif
  
  endwhile
  
  tagE=p;
  
endfunction

%[xml, tagC, tagE]=parseXMLTag(txt)
function [xml, tagC, tagE]=parseXMLTag(txt)
  
  xml=struct;
  tagC=false;
  
  %tag name
  [S,E, TE]=regexp(txt, '^<(([a-zA-Z0-9]+):)*([a-zA-Z0-9]+)\s*');
  if( size(TE{1},1) > 1 )
    xml.__prefix=txt(TE{1}(2,1):TE{1}(2,2));
    xml.__name=txt(TE{1}(3,1):TE{1}(3,2));
  else
    xml.__name=txt(TE{1}(1,1):TE{1}(1,2));
  endif


  p=1+E;
  while( p < length(txt) )

    [S,E, TE]=regexp(txt(p:end), ['^([a-zA-Z0-9]+)\s*=\s*("([^"]+)"|' "'" '([^' "'" ']+)' "'" ')\s*']);
    if( any(S) )
      xml=setfield(xml, [ '_' txt(p-1+[TE{1}(1,1):TE{1}(1,2)]) ], txt(p-1+[TE{1}(3,1):TE{1}(3,2)]));
      p=p+E;
      continue
    endif
    
    [S,E]=regexp(txt(p:end), '^/>');
    if( any(S) )
      tagC=true;
      tagE=p+E;
      return
    endif
    
    [S,E]=regexp(txt(p:end), '^>');
    if( any(S) )
      tagC=false;
      tagE=p+E;
      return
    endif
    
    %syntax error: skip one character
    p=p+1;
    
  endwhile

  tagE=p;
  
endfunction
