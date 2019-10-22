clearvars; close all; clc; 
% Read in the Marmousi ii p-wave velocity model.
% Downloaded from here: http://www.agl.uh.edu/downloads/downloads.htm
%

%-----------------------------------------------------------
%   Keith Roberts   : 2019 --
%   Email           : krober@usp.br
%   Last updated    : 10/20/2019
%-----------------------------------------------------------
%

% ensure path is set correctly
libpath 

fname = 'MODEL_P-WAVE_VELOCITY_1.25m.segy'; 

gdat = geodata('segy',fname,'gridspace',1.25) ; 

%plot(gdat) % visualize p-wave velocity model

fh = edgefx('wl',500,'geodata',gdat,'min_el',400,'max_el',5e3,'g',0.50); 

%plot(fh); % visualize mesh size function

mshopts = meshgen('geodata',gdat,'edgefx',fh,'plot_on',1);





