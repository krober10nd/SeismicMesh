clearvars; close all; clc;

% Read in simple synthetic 3D P-wave model to build surface mesh .
%
%-----------------------------------------------------------
%   Keith Roberts   : 2019 --
%   Email           : krober@usp.br
%   Last updated    : 11/21/2019
%-----------------------------------------------------------
%

% ensure path is set correctly
libpath

%% specify mesh design
FNAME  = 'SeismicUnitCubeVp.nc'; % file containing velocity model
OFNAME = 'SeismicUnitCubeVp'; % output file name
MIN_EL = 100 ; % minimum element size (meters)
MAX_EL = 5e3 ;% maximum element size (meters)
WL     = 100 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 20 ; % maximum shot record frequency (hz)
GRADE  = .15 ; % expansion rate of element size
ITMAX  = 30 ;
%% build sizing function
gdat = geodata('velocity_model',FNAME) ;

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,...
    'min_el',MIN_EL,'max_el',MAX_EL,'g',GRADE);
      

% control points for a cubic surface
cp{1,1}=[0,0,0];cp{1,2}=[0,0.33,0];cp{1,3}=[0,0.66,0];cp{1,4}=[0,1.0,0];

cp{2,1}=[0.33 0 0];cp{2,2}=[0.33 0.33 0]; cp{2,3}=[0.33 0.66 0]; cp{2,4}=[0.33 1 0];

% original
%cp{3,1}=[0.66 0 0];cp{3,2}=[0.66 0.33 0]; cp{3,3}=[0.66 0.66 0]; cp{3,4}=[0.66 1 0]; 
% modified
cp{3,1}=[0.66 0 0];cp{3,2}=[0.66 0.33 0]; cp{3,3}=[0.66 0.66 3]; cp{3,4}=[0.66 1 0]; 


cp{4,1}=[1 0 0];cp{4,2}=[1 0.33 0];cp{4,3}=[1 0.66 0];cp{4,4}=[1 1 0];


w=cell(4,4); % weights
for i=1:4
    for j=1:4
        cp{i,j}(1)=cp{i,j}(1)*4900;
        cp{i,j}(2)=cp{i,j}(2)*4900;
        cp{i,j}(3)=cp{i,j}(3)*4900;
        w{i,j}(1)=1;
        w{i,j}(2)=1;
        w{i,j}(3)=1;
    end
end

fh = @(p) ef.F(p) ;

[p,t]=meshsurface(cp,w,fh,MIN_EL,gdat.bbox',ITMAX);

%%
pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids

tria = triangulation(t,p); 

%% for each face, place a node 
for i =  1 : tria.size(1)
    fN=tria.faceNormal(i);
    pause
end
