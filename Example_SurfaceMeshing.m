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
MIN_EL = 60 ; % minimum element size (meters)
MAX_EL = 5e3 ;% maximum element size (meters)
WL     = 20 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 20 ; % maximum shot record frequency (hz)
GRADE  = 2 ; % expansion rate of element size
ITMAX  = 30 ;
%% build sizing function
gdat = geodata('velocity_model',FNAME) ;

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,...
    'min_el',MIN_EL,'max_el',MAX_EL,'g',GRADE);
      
% control points for a quadratic surface 
cp{1,1}=[0,0,0];cp{1,2}=[0,0.5,0];cp{1,3}=[0,1,0];
cp{2,1}=[0.5,0,0.5];
cp{2,2}=[0.2,0.0,2.0]; % should be 0.5 0.5 0.5 for a linear surface
cp{2,3}=[0.5,1,0.5];
cp{3,1}=[1,0,1];cp{3,2}=[1,0.5,1];cp{3,3}=[1,1,1];
w=cell(3,3); % weights
for i=1:3
    for j=1:3
        cp{i,j}(1)=cp{i,j}(1)*4900;
        cp{i,j}(2)=cp{i,j}(2)*4900;
        cp{i,j}(3)=cp{i,j}(3)*4900;
        
        w{i,j}(1)=1;
        w{i,j}(2)=1;
        w{i,j}(3)=1;

    end
end

w{2,2}(2)=1;


fh = @(p) ef.F(p) ;

[p,t]=meshsurface(cp,w,fh,MIN_EL,gdat.bbox',ITMAX);

