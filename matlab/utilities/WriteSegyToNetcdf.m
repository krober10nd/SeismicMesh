clearvars; close all; clc; 
%%
OFNAME = 'SeismicUnitCubeVp_V2.nc';
data=ncread('SeismicUnitCubeVp.nc','vp'); 
dims = size(data);
[nx,ny,nz]=size(data);

%% create attributes of file 


nccreate(OFNAME,'dim','Dimensions',{'dim',1});
ncwriteatt(OFNAME, 'dim', 'long_name', 'dimension of problem');
ncwriteatt(OFNAME,'dim','units','dimension'); 



nccreate(OFNAME,'gridspace','Dimensions',{'gridspace',1});
ncwriteatt(OFNAME, 'gridspace', 'long_name', 'grid spacing in meters');
ncwriteatt(OFNAME,'gridspace','units','meters'); 



nccreate(OFNAME,'nx','Dimensions',{'nx',1});
ncwriteatt(OFNAME, 'nx', 'long_name', 'number of coordinates in x-direction');
ncwriteatt(OFNAME,'nx','units','gridpoints'); 


nccreate(OFNAME,'ny','Dimensions',{'ny',1});
ncwriteatt(OFNAME, 'ny', 'long_name', 'number of coordinates in y-direction');
ncwriteatt(OFNAME,'ny','units','gridpoints'); 


nccreate(OFNAME,'nz','Dimensions',{'nz',1});
ncwriteatt(OFNAME, 'nz', 'long_name', 'number of coordinates in z-direction');
ncwriteatt(OFNAME,'nz','units','gridpoints'); 


nccreate(OFNAME,'x0y0z0','Dimensions',{'x0y0z0',length(dims)});
ncwriteatt(OFNAME, 'x0y0z0', 'long_name', 'position of bottom left front corner');
ncwriteatt(OFNAME,'x0y0z0','units','cartesian coordinates'); 


nccreate(OFNAME,'vp','Dimensions',{'x',50,'y',50,'z',50});
ncwriteatt(OFNAME, 'vp', 'long_name', 'P-wave velocity');
ncwriteatt(OFNAME,'vp','units','m/s')

%% write data to file 

ncwrite(OFNAME,'dim',length(dims));

ncwrite(OFNAME,'gridspace',100);

ncwrite(OFNAME,'x0y0z0',[0,0,0]);

ncwrite(OFNAME,'nx',dims(1));

ncwrite(OFNAME,'ny',dims(2));

ncwrite(OFNAME,'nz',dims(3));


ncwrite(OFNAME,'vp',data);



ncdisp(OFNAME);



