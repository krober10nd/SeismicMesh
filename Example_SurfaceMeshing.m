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
MIN_EL = 40 ; % minimum element size (meters)
MAX_EL = 5e3 ;% maximum element size (meters)
WL     = 20 ;% number of nodes per wavelength of wave with FREQ (hz)
FREQ   = 20 ; % maximum shot record frequency (hz)
GRADE  = 2 ; % expansion rate of element size
ITMAX  = 10 ;
%% build sizing function
gdat = geodata('velocity_model',FNAME) ;

ef = edgefx('geodata',gdat,...
    'wl',WL,'f',FREQ,...
    'min_el',MIN_EL,'max_el',MAX_EL,'g',GRADE);
      
% control points for a linear surface 
cp{1,1} = [0 0 0];    cp{1,2} = [0 4900 4900];
cp{2,1} = [4900 0 0]; cp{2,2} = [4900 4900 4900];
    
fh = @(p) ef.F(p) ;

[p,t]=meshsurface(cp,fh,MIN_EL,gdat.bbox',ITMAX);

