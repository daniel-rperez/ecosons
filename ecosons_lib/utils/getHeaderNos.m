%funcion auxiliar para extraer los numeros de las columnas de interes en un CSV
% c=getHeaderNos(h,...)
% c: vector con los números de las columnas seleccionadas
% h: cadena con los nombres de las columnas separados por comas
% ...: cadenas con los nombres de las columnas que se quieren seleccionar
%%Version: 2026/08/22
function c=getHeaderNos(h,varargin)
	if( ischar(h) )
		h=strsplit(h, ',');
	endif
	
	c=[];
	for n=1:length(varargin)
		m=find(strcmp(h,varargin{n}));
		if( any(m) )
			c(n)=m;
		else
			c(n)=0;
		endif
	endfor
endfunction
