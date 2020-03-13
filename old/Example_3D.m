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
FNAME  = 'SeismicUnitCubeVp.nc'; % file containing velocity model
OFNAME = 'SeismicUnitCubeVp'; % output file name
IN_EL = 50 ; % minimum element size (meters)
MAX_EL = 5e3 ;% maximum element size (meters)
WL     = 50 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 3 ; % maximum shot record frequency (hz)
GRADE  = 0.35 ; % expansion rate of element size
SLP    = 50;
%% build sizing function
gdat = geodata('velocity_model',FNAME) ;

plot(gdat);

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,'slp',SLP,...
    'min_el',MIN_EL,'max_el',MAX_EL,'g',GRADE);

plot(ef);



%% mesh generation step

drectangle3D=@(p,x1,x2,y1,y2,z1,z2)-min(min(min(min(min(-z1+p(:,3),z2-...
    p(:,3)),-y1+p(:,2)),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

fd = @(p) drectangle3D(p,...
    gdat.bbox(1,1),gdat.bbox(1,2),...
    gdat.bbox(2,1),gdat.bbox(2,2),...
    gdat.bbox(3,1),gdat.bbox(3,2));

fh = @(p) ef.F(p);

P_FIX=[]; % INSERT FIXED POINTS HERE 
IT_MAX=100; % DEFAULT NUMBER OF MESH GENERATION ITERATIONS 1000

tic
[P,T,STATS] = distmesh( fd, fh, MIN_EL, gdat.bbox', [] , [], IT_MAX );
toc

%% write the mesh to disk (0,0) is top left not bottom left corner. 
% flip it upside down 
P(:,2)=P(:,2)*-1;

%% visualize result
figure; simpplot(P,T) ; 

%% write mesh to a msh format
gmsh_mesh2d_write ( 'output.msh', 2, length(P), P', ...
  3, length(T), T' ) ;

vtkwrite('volume_wsurface.vtk','polydata','tetrahedron',P(:,1),P(:,2),P(:,3),T)
