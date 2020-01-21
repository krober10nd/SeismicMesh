
clearvars; close all; clc;
% Write the mesh file's points and vertices to a gmsh file. 

%-----------------------------------------------------------
%   Keith Roberts   : 2019 --
%   Email           : krober@usp.br
%   Last updated    : 1/19/2020
%-----------------------------------------------------------
%

nameOfMsh = 'Testing.msh' ; 
load MM_Subdomain_SMALL
NP = length(P3); % number of nodes
NT = length(T3); % number of cells

keep = zeros(NT,1); 

% NB: we use version 2.2 of the gmsh format. 
% the header looks like this...
% $MeshFormat
% 2.2 0 8
% $EndMeshFormat

fileID = fopen(nameOfMsh,'w');
fprintf(fileID,'%s\n','$MeshFormat');
fprintf(fileID,'%s\n','2.2 0 8');
fprintf(fileID,'%s\n','$EndMeshFormat');

fprintf(fileID,'%s\n','$Nodes'); 
fprintf(fileID,'%i\n',NP); 

disp('Writing nodal table...')
for i = 1 : NP
   nodal_table(i,1) = i ; 
   nodal_table(i,2) = P3(i,1)./1e3 ; % convert to km
   nodal_table(i,3) = P3(i,2)./1e3 ; % convert to km
   nodal_table(i,4) = 0 ; 
   
   fprintf(fileID,'%i %f %f %i\n',nodal_table(i,:)); 
   
end
fprintf('%s\n','$EndNodes'); 

disp('Writing element table...'); 
fprintf(fileID,'%s\n','$Elements'); 
fprintf(fileID,'%i\n',NT); 
for i = 1 : NT
    if keep(i)
        element_table(i,1) = i; % element number
        element_table(i,2) = 2; % type of element 
        element_table(i,3) = 2; % number of tags 
        element_table(i,4) = 4; % tag of physical thing
        element_table(i,5) = 1; % tag of element thing
        
        element_table(i,6) = T3(i,1) ;
        element_table(i,7) = T3(i,2) ;
        element_table(i,8) = T3(i,3) ;
        
        fprintf(fileID,'%i %i %i %i %i %i %i %i\n',element_table(i,:)); 

    else
        element_table(i,1) = i; % element number
        element_table(i,2) = 2; % type of element
        element_table(i,3) = 2; % number of tags
        element_table(i,4) = 3; % tag of physical thing
        element_table(i,5) = 2; % tag of element thing
        
        element_table(i,6) = T3(i,1) ;
        element_table(i,7) = T3(i,2) ;
        element_table(i,8) = T3(i,3) ;
        
        fprintf(fileID,'%i %i %i %i %i %i %i %i\n',element_table(i,:));
        
    end
end
fprintf('%s\n','$EndElements'); 

fclose(fileID);

disp('Finished writing file...'); 

figure; simpplot(P3,T3); 
