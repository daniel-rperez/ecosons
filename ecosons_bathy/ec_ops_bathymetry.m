%[ntr, utmCoords, xCoord, yCoord, znCoord, depth]=ec_ops_bathymetry
%
%
function [ntr, utmCoords, xCoord, yCoord, znCoord, depth, bTime]=ec_ops_bathymetry
 global BATHYMETRY

 %null values
 ntr=[];
 utmCoords=false;
 xCoord=[];
 yCoord=[];
 znCoord=[];
 depth=[];
 bTime=[];

 if( length(BATHYMETRY)==0 )
  return;
 endif

 nn=1;
 lat=[];
 lon=[];
  
 for sds=1:length(BATHYMETRY)
  ll=length(BATHYMETRY{sds}.time);

  ntr(nn:nn+ll-1)=sds;

  lat(nn:nn+ll-1)=BATHYMETRY{sds}.latitude;
  lon(nn:nn+ll-1)=BATHYMETRY{sds}.longitude;

  depth(nn:nn+ll-1)=BATHYMETRY{sds}.depth;

  bTime(nn:nn+ll-1)=BATHYMETRY{sds}.time;
  
  nn=nn+ll;
 endfor
  
 %convert to UTM? (default no)
 ask=gminput('Use UTM coordinates? (y/N) ', 's');
 if( length(ask)~=0 && ( ask(1)=='y' || ask(1)=='Y' ) )
  
  utmCoords=true;
  
  znCoord=round((lon(1)+183)/6); %use first point zone
  
  [xCoord, yCoord, znCoord]=latlon2utmxy(znCoord, lat,lon);
  
 else
  utmCoords=false;
  xCoord=lon;
  yCoord=lat;
  znCoord=NaN;
 endif

endfunction
