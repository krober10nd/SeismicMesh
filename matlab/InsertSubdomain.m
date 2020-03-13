
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
OFNAME = 'Marmousi_subse_modified_1_immersed'; % output fi
MIN_EL = 10 ; % minimum element size (meters)
MAX_EL = 1e3 ;% maximum element size (meters)
WL     = 20 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 5 ; % maximum shot record frequency (hz)
SLP    = 0 ; % element size (meters) near maximum gradient in P-wavepseed.
GRADE  = 0.10 ; % expansion rate of mesh resolution (in decimal percent).
CR     = 0.3 ; % desired Courant number desired DT will be enforced for.
DT     = 2e-3; % desired time step (seconds).
%% build sizing function
gdat = geodata('velocity_model',FNAME) ;

% plot(gdat) % visualize p-wave velocity model

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,...
    'dt',DT,'cr',CR,...
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

IT_MAX=100; % NUMBER OF MESH GENERATION ITERATIONS 

load Subdomain2
% create 3 subdomains spread across the top of the domain by shifting the
% first one

num_subdomains = 5 ;
PshiftZ150 = [P2(:,1), P2(:,2) + 150];
PshiftZ300 = [P2(:,1), P2(:,2) + 300];
PshiftZ450 = [P2(:,1), P2(:,2) + 450];
PshiftZ600 = [P2(:,1), P2(:,2) + 600];
PshiftZ750 = [P2(:,1), P2(:,2) + 750];

PFIX  =  [PshiftZ50 ; PshiftZ300; PshiftZ450; PshiftZ600; PshiftZ750];
T2 = TFIX; 
for i = 1 : num_subdomains
    EGFIX = [EE2; EE2 + i*length(P2);];
    TFIX  = [TFIX ; T2 + length(P2)*i;];
end
[PFIX,TFIX]=fixmesh(PFIX,TFIX); 
BNDE = extdom_edges2(TFIX,PFIX);
POLY_FIX = cell2mat(extdom_polygon(BNDE,PFIX,-1)');

[P,T,STATS] = distmesh( fd, fh, MIN_EL, gdat.bbox', PFIX, EGFIX, POLY_FIX, IT_MAX);
%[P,T,STATS] = distmesh( fd, fh, MIN_EL, gdat.bbox', [], [], IT_MAX );

%% Mesh improvement 

[P3,T3]=direct_smoother_lur(P,T,PFIX,1);

%% Visualization save MM_Subdomain P3 T3 P2

figure; simpplot(P3,T3)

save MM_Subdomain_SMALL P3 T3 TFIX PFIX 
