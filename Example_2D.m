clearvars; close all;% clc;
% Read in the Marmousi ii p-wave velocity model.
% Downloaded from here: http://www.agl.uh.edu/downloads/downloads.htm
%

%-----------------------------------------------------------
%   Keith Roberts   : 2019 --
%   Email           : krober@usp.br
%   Last updated    : 11/3/2019
%-----------------------------------------------------------
%

% ensure path is set correctly
libpath
%% specify mesh design
FNAME = 'MODEL_P-WAVE_VELOCITY_1.25m.segy'; % SegY file containing velocity model
GRIDSPACE = 1.25 ; % grid spacing (meters) p-wave model.
OFNAME = 'OUTPUT.msh'; % output fi
MIN_EL = 10 ; % minimum element size (meters)
MAX_EL = 5e3 ;% maximum element size (meters)
WL     = 0 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 0 ; % maximum shot record frequency (hz)
SLP    = 50 ; % element size (meters) near maximum gradient in P-wavepseed.
GRADE  = 0.50 ; % expansion rate of mesh resolution (in decimal percent).
CR     = 0.1 ; % desired Courant number desired DT will be enforced for.
DT     = 1e-3; % desired time step (seconds).
%% build sizing function
gdat = geodata('segy',FNAME,'gridspace',GRIDSPACE) ;

plot(gdat) % visualize p-wave velocity model

ef = edgefx('geodata',gdat,...
    'slp',SLP,...
    'wl',WL,'f',FREQ,...
    'dt',DT,'cr',CR,...
    'min_el',MIN_EL,'max_el',MAX_EL,...
    'g',GRADE);

 plot(ef); % visualize mesh size function

%% mesh generation step
drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

fd = @(p) max( drectangle(p,...
     gdat.bbox(1,1),gdat.bbox(1,2),gdat.bbox(2,1),gdat.bbox(2,2)),...
     -(sqrt(sum(p.^2,2))-0.5) );
 
fh = @(p) ef.F(p); 

P_FIX=[]; % INSERT FIXED POINTS HERE 
IT_MAX=100; % NUMBER OF MESH GENERATION ITERATIONS 

[ P, T, COUNT ] = distmeshnd( fd, fh, MIN_EL, gdat.bbox', P_FIX, IT_MAX ) ;


%% Put the velocity onto the nodes and write this to disk 
[vp_onNodes] = gdat.interpVp(P(:,1:2)); 
% write out the vp file 
fid = fopen('MarmousII_VpOnNodes.txt','w'); 
for i = 1 : length(vp_onNodes) 
   fprintf(fid,'%f\n',vp_onNodes(i));  
end
fclose(fid); 
%% write the mesh to disk (0,0) is top left not bottom left corner. 
% flip it upside down 
P(:,2)=P(:,2)*-1;

%% visualize result
figure; patch( 'vertices', P, 'faces', T, 'facecolor', [.90 .90 .90] )

axis equal
axis tight

figure; trisurf(T,P(:,1),P(:,2),vp_onNodes); shading interp; 
axis equal; axis tight; view(2); 


%% write mesh to a msh format
gmsh_mesh2d_write ( 'output.msh', 2, length(P), P', ...
  3, length(T), T' ) ;


% write mesh to a vtu format
P(:,3) = zeros*P(:,2); % add dummy 3rd dimension
exportTriangulation2VTK('output',P,T) ;


