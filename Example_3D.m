clearvars; close all; clc;
% Read in the UnitCube.segy a simple synthetic 3D P-wave velocity model.
%
%-----------------------------------------------------------
%   Keith Roberts   : 2019 --
%   Email           : krober@usp.br
%   Last updated    : 10/28/2019
%-----------------------------------------------------------
%

% ensure path is set correctly
libpath

%% specify mesh design
FNAME = 'UnitCubeVp.nc'; % file containing velocity model
GRIDSPACE = 50 ; % grid spacing (meters) p-wave model.
OFNAME = 'UnitCube'; % output file name
MIN_EL = 25 ; % minimum element size (meters)
MAX_EL = 5e3 ;% maximum element size (meters)
WL     = 500 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 15 ; % maximum shot record frequency (hz)
GRADE  = 0.35 ; % expansion rate of element size
DIM    = 3 ; % 3D problem!
%% build sizing function
gdat = geodata('segy',FNAME,'dim',DIM,'gridspace',GRIDSPACE) ;

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,...
    'min_el',MIN_EL,'max_el',MAX_EL);

%% mesh generation step

drectangle3D=@(p,x1,x2,y1,y2,z1,z2)-min(min(min(min(min(-z1+p(:,3),z2-p(:,3)),-y1+p(:,2)),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

fd = @(p) drectangle3D(p,...
    gdat.bbox(1,1),gdat.bbox(1,2),...
    gdat.bbox(2,1),gdat.bbox(2,2),...
    gdat.bbox(3,1),gdat.bbox(3,2));

fh = @(p) ef.F(p);

P_FIX=[]; % INSERT FIXED POINTS HERE 
IT_MAX=100; % DEFAULT NUMBER OF MESH GENERATION ITERATIONS 1000

[ P, T, COUNT ] = distmeshnd( fd, fh, MIN_EL, gdat.bbox', P_FIX, IT_MAX ) ;
P(:,2)=P(:,2)*-1;

F = [T(:,[1:3]); T(:,[1,2,4]); T(:,[2,3,4]); T(:,[3,1,4])];

figure; patch( 'vertices', P, 'faces', F, 'facecolor', [.9, .9, .9] ), view(3)

axis equal
