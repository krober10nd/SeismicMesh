clearvars; close all; clc;
% Read in the Marmousi ii p-wave velocity model.
% Downloaded from here: http://www.agl.uh.edu/downloads/downloads.htm
%

%-----------------------------------------------------------
%   Keith Roberts   : 2019 --
%   Email           : krober@usp.br
%   Last updated    : 10/20/2019
%-----------------------------------------------------------
%

% ensure path is set correctly
libpath
%%
MIN_EL = 10 ; 
MAX_EL = 5e3 ;
WL     = 10 ; 
GRADE  = 0.90 ; 
GRIDSPACE = 1.25 ; 
FNAME = 'MODEL_P-WAVE_VELOCITY_1.25m.segy'; 
%%
gdat = geodata('segy',FNAME,'gridspace',GRIDSPACE) ;

%plot(gdat) % visualize p-wave velocity model

ef = edgefx('wl',WL,'geodata',gdat,'min_el',MIN_EL,'max_el',MAX_EL,'g',GRADE);

%plot(fh); % visualize mesh size function

%% mesh generation step
drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

fd = @(p) max( drectangle(p,...
     gdat.bbox(1,1),gdat.bbox(1,2),gdat.bbox(2,1),gdat.bbox(2,2)),...
     -(sqrt(sum(p.^2,2))-0.5) );
 
fh = @(p) ef.F(p); 

P_FIX=[]; 
E_FIX=[]; 
IT_MAX=100; % DEFAULT 1000
FID=1;
FIT=[];

[ P, T, STAT ] = distmesh( fd, fh, MIN_EL, gdat.bbox', P_FIX, E_FIX, IT_MAX, FID, FIT ) ;


patch( 'vertices', P, 'faces', T, 'facecolor', [.9, .9, .9] )




