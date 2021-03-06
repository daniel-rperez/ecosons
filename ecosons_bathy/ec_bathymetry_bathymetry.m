%[err,err_desc]=ec_bathymetry_bottom()
%creates a void bathymetry from selected transects
%no arguments required
%returns error code and description
%CALLS: 
function [err, err_desc]=ec_bathymetry_bathymetry
 global SONAR_DATA
 global SONAR_DATA_SELECTION
 global BATHYMETRY
 err=0;
 err_desc='';

 %check data availability
 if( ~iscell(SONAR_DATA) || length(SONAR_DATA)==0 )
  err=1;
  err_desc='Error: no echograms defined';
  return;
 endif

 %check if echosounder provides a bathymetry
 sounder_has_depth=0;
 for sds=SONAR_DATA_SELECTION
  sounder_has_depth=sounder_has_depth || isfield(SONAR_DATA{sds}.Q, 'depth');
 endfor
 if( sounder_has_depth )
  ask=gminput('Use echosounder bathymetry? [y/N] ', 's');
  if( length(ask)==0 || ~( ask(1)=='y' || ask(1)=='Y') )
   sounder_has_depth=0;
  endif
 endif

 %ask for bathymetry cutoff
 depth0=gminput('Minimum depth cutoff (m, default none): ');
 if( length(depth0)~=1 )
  depth0=0;
 endif

 %bathymetry structure cells (one per transect)
 BATHYMETRY={};
 nt=0;

 for sds=SONAR_DATA_SELECTION
  dta=SONAR_DATA{sds};
 
  nn=0;
  lat=[];
  lon=[];
  tme=[];
  dep=[];
  
  for n=1:length(dta.G)
 
   if( sounder_has_depth )
    dd=dta.Q(n).depth; %use echosounder depth
   else
    dd=0.5*dta.Q(n).soundVelocity*(dta.R(n)-2)*dta.Q(n).sampleInterval;
   endif
   
   if( dta.G(n).time<0 || dd<depth0 )
    continue;
   endif

   nn=nn+1;
   lat(nn)=dta.G(n).latitude;
   lon(nn)=dta.G(n).longitude;
   tme(nn)=dta.G(n).time;
   dep(nn)=dd;
  
  endfor

 %fill bathymetry object cells
 if( length(tme)>0 )
  nt=nt+1;
  BATHYMETRY{nt}.name=dta.name;
  BATHYMETRY{nt}.latitude=lat;
  BATHYMETRY{nt}.longitude=lon;
  BATHYMETRY{nt}.time=tme;
  BATHYMETRY{nt}.depth=dep;
 endif

 endfor

endfunction


