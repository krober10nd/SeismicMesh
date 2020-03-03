clearvars; close all; clc;

libpath
%% specify mesh design
FNAME  = 'Full_Marmousi_Initial_Guess.nc'; % file containing velocity model
OFNAME = 'Marmousi_full_initial_DOMAIN_EXT.msh'; % output fi
WL     = 20 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 20 ; % maximum shot record frequency (hz)
%% build sizing function
gdat = geodata('velocity_model',FNAME,'expand',500) ;

%plot(gdat) % visualize p-wave velocity model

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,'max_sponge_res',5e3);

%plot(ef); % visualize mesh size function

%% mesh generation step
drectangle=@(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

fd = @(p) max( drectangle(p,...
     gdat.bbox(1,1),gdat.bbox(1,2),gdat.bbox(2,1),gdat.bbox(2,2)),...
     -(sqrt(sum(p.^2,2))-0.5) );
 
fh = @(p) ef.F(p); 

P_FIX=[]; % INSERT FIXED POINTS HERE 
IT_MAX=100; % NUMBER OF MESH GENERATION ITERATIONS 

[P,T,STATS] = distmesh( fd, fh, 25, gdat.bbox', [], [], [], IT_MAX );

%% Mesh improvement 

[P2,T2]=direct_smoother_lur(P,T,[],1);
figure; triplot(T2,P2(:,1),P2(:,2))

P2=protate(P2,(90*pi)/180);
P2(:,1)=P2(:,1)./1e3; % in kms
P2(:,2)=P2(:,2)./1e3; % in kms

figure; triplot(T2,P2(:,1),P2(:,2))

%% write mesh to a msh format
gmsh_mesh2d_write ( OFNAME, 2, length(P2), P2', ...
  3, length(T2), T2' ) ;