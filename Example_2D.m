clearvars; close all; clc;
% Read in a subsection of the Marmousi ii p-wave velocity model.
% File is included in the data directory. 

%-----------------------------------------------------------
%   Keith Roberts   : 2019 --
%   Email           : krober@usp.br
%   Last updated    : 11/19/2019
%-----------------------------------------------------------
%

% ensure path is set correctly
libpath
%% specify mesh design
FNAME  = 'Marmousi_subset_modified_1.nc'; % file containing velocity model
OFNAME = 'Marmousi_subse_modified_1'; % output fi
MIN_EL = 10 ; % minimum element size (meters)
MAX_EL = 1e3 ;% maximum element size (meters)
WL     = 20 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 20 ; % maximum shot record frequency (hz)
SLP    = 10 ; % element size (meters) near maximum gradient in P-wavepseed.
GRADE  = 0.10 ; % expansion rate of mesh resolution (in decimal percent).
CR     = 0.1 ; % desired Courant number desired DT will be enforced for.
DT     = 1e-3; % desired time step (seconds).
%% build sizing function
gdat = geodata('velocity_model',FNAME) ;

% plot(gdat) % visualize p-wave velocity model

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,...
    'dt',DT,'cr',CR,...
    'slp',SLP,...
    'g',GRADE,...
    'min_el',MIN_EL,...
    'max_el',MAX_EL);

plot(ef); % visualize mesh size function

%% mesh generation step
drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

fd = @(p) max( drectangle(p,...
     gdat.bbox(1,1),gdat.bbox(1,2),gdat.bbox(2,1),gdat.bbox(2,2)),...
     -(sqrt(sum(p.^2,2))-0.5) );
 
fh = @(p) ef.F(p); 

P_FIX=[]; % INSERT FIXED POINTS HERE 
IT_MAX=50; % NUMBER OF MESH GENERATION ITERATIONS 

[P,T,STATS] = distmesh( fd, fh, MIN_EL, gdat.bbox', [], [], IT_MAX );

%% Mesh improvement 

[P2,T2]=direct_smoother_lur(P,T,[],1);

vp_on_nodes = InterpVP(gdat,P2); 

P2=protate(P2,(90*pi)/180);

figure; trisurf(T2,P2(:,1),P2(:,2),vp_on_nodes); shading interp; view(2); 

fid = fopen('VpOnNodes.txt','w'); 
for i = 1 : length(vp_on_nodes)
   fprintf(fid,'%f\n',vp_on_nodes(i)); 
end
fclose(fid); 
%% Visualization 
figure; simpplot(P2,T2)

PlotMeshResolution(P2,T2); 
%% write mesh to a msh format
gmsh_mesh2d_write ( OFNAME, 2, length(P2), P2', ...
  3, length(T2), T2' ) ;

% write mesh to a vtu format
P2(:,3) = zeros*P2(:,2); % add dummy 3rd dimension
exportTriangulation2VTK(OFNAME,P2,T2) ;
 

