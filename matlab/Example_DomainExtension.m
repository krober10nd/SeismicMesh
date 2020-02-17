
clearvars; close all; clc;
% Read in a subsection of the Marmousi ii p-wave velocity model.
% File is included in the data directory. 
% --Build a domain extension around the border of the domain

%-----------------------------------------------------------
%   Keith Roberts   : 2019 --
%   Email           : krober@usp.br
%   Last updated    : 1/19/2020
%-----------------------------------------------------------
%

% ensure path is set correctly
libpath
%% specify mesh design
FNAME  = 'Marmousi_subset_modified_1.nc'; % file containing velocity model
OFNAME = 'Marmousi_subse_modified_1_immersed'; % output fi
MIN_EL = 10 ; % minimum element size (meters)
MAX_EL = 500 ;% maximum element size (meters)
WL     = 20 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 5 ; % maximum shot record frequency (hz)
SLP    = 0 ; % element size (meters) near maximum gradient in P-wavepseed.
GRADE  = 3.0 ; % expansion rate of mesh resolution (in decimal percent).
CR     = 0.3 ; % desired Courant number desired DT will be enforced for.
DT     = 2e-3; % desired time step (seconds).
max_sponge_res = 500; % maximum element size in domain extension (meters)
expand = 2e3; % amount in x and z directions to expand domain
%% build sizing function
gdat = geodata('velocity_model',FNAME,'expand',2e3) ;

% plot(gdat) % visualize p-wave velocity model

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,...
    'dt',DT,'cr',CR,...
    'min_el',MIN_EL,...
    'max_el',MAX_EL,...
    'max_sponge_res',max_sponge_res);

plot(ef); % visualize mesh size function

%% mesh generation step
drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

fd = @(p) max( drectangle(p,...
     gdat.bbox(1,1),gdat.bbox(1,2),gdat.bbox(2,1),gdat.bbox(2,2)),...
     -(sqrt(sum(p.^2,2))-0.5) );
 
fh = @(p) ef.F(p); 

IT_MAX=200; % NUMBER OF MESH GENERATION ITERATIONS 

[P,T,STATS] = distmesh( fd, fh, MIN_EL, gdat.bbox', [], [], [], IT_MAX );

%% Mesh improvement 

[P3,T3]=direct_smoother_lur(P,T,[],1);

P3=protate(P3,(90*pi)/180);

%% Visualization save MM_Subdomain P3 T3 P2

figure; simpplot(P3,T3)

save MM_Subdomain_SMALL P3 T3
