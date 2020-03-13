clearvars; close all; clc;
%% Make a small strip along the top of the domain to be constrained
% in the larger mesh.
libpath

FNAME  = 'Marmousi_subset_modified_1.nc'; % file containing velocity model
OFNAME = 'Marmousi_subse_modified_1'; % output fi
MIN_EL = 10 ; % minimum element size (meters)
MAX_EL = 1e3 ;% maximum element size (meters)
WL     = 20 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 5 ; % maximum shot record frequency (hz)
GRADE  = 0.10 ; % expansion rate of mesh resolution (in decimal percent).
CR     = 0.3 ; % desired Courant number desired DT will be enforced for.
DT     = 2e-3; % desired time step (seconds).
%% build sizing function
gdat = geodata('velocity_model',FNAME) ;

bbox = [ 50  100
         50  100];
    
bbox = [50 100 % z
        50 100]; % x
    
%bbox=[           0        2000
%                 0        3000];            
% L = 3000; 
% W = 2000; 

% ytop = 0.04*2000 + 25; 
% ybot = ytop - 0.025*W;
% 
% xleft = 0.05*W;
% xright = xleft + L-(0.05*L); 
% 
% sub_bbox = [ybot ytop 
%             xleft xright]; 
%gdat.bbox = sub_bbox; 

gdat.bbox = bbox; 

% plot(gdat) % visualize p-wave velocity model

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,...
    'dt',DT,'cr',CR,...
    'g',GRADE,...
    'min_el',MIN_EL,...
    'max_el',MAX_EL);

%plot(ef); % visualize mesh size function


drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

fd = @(p) max( drectangle(p,...
     gdat.bbox(1,1),gdat.bbox(1,2),gdat.bbox(2,1),gdat.bbox(2,2)),...
     -(sqrt(sum(p.^2,2))-0.5) );
 
fh = @(p) ef.F(p); 

P_FIX=[]; % INSERT FIXED POINTS HERE 
IT_MAX=500; % NUMBER OF MESH GENERATION ITERATIONS 

[P,T,STATS] = distmesh( fd, fh, MIN_EL, gdat.bbox', [], [], IT_MAX );


[P2,T2]=direct_smoother_lur(P,T,[],1);

%% Turn mesh into point and edge constraints 
EE2 = [T2(:,[1,2]); T2(:,[1,3]); T2(:,[2,3])];           % Interior bars duplicated
EE2 = unique(sort(EE2,2),'rows');                    % Bars as node pairs

EE2 = renumberEdges(EE2) ; 

save Subdomain2 P2 T2 EE2

