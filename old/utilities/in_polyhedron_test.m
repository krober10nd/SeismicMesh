%% Tutorial and tests of IN_POLYHEDRON function
% *By Jarek Tuszynski* (jaroslaw.w.tuszynski@leidos.com)
%
% IN_POLYHEDRON tests if points are inside a 3D triangulated surface 
% (faces/vertices) or volume (tetrahedrals/vertices). There are NO
% assumptions about orientation of the face normals.
%
% IN = INPOLYHEDRON(X,POINTS) tests if the query points (POINTS) are inside
%   the surface/polyhedron defined by X. X can be a structure with fields 
%   'vertices' and 'faces' or an object of MATLAB triangulation class. In
%   case of triangulation class object we will only use the outside
%   boundary. POINTS is an N-by-3 set of XYZ coordinates. IN is an N-by-1 
%   logical vector which will be TRUE for each query point inside the surface.
%
% INPOLYHEDRON(FACES,VERTICES,POINTS) takes faces/vertices separately, 
%   rather than in an FV structure.
%
%% Algorithm
% For each point do:
%
% # shoot a random ray out of the query point in a random direction
% # for each face solve:
% $\left[\begin{array}{ccc} -d_{x} & v1_{x}-v0_{x} & v2_{x}-v0_{x} \\ -d_{y} & v1_{y}-v0_{y} & v2_{y}-v0_{y} \\ -d_{z} & v1_{z}-v0_{z} & v2_{z}-v0_{z} \end{array}\right]\*\left[\begin{array}{c} t \\ u \\ v \end{array} \right]=\left[\begin{array}{c} o_{x}-v0_{x} \\ o_{y}-v0_{y} \\ o_{z}-v0_{z} \end{array}\right]$
% for $\left[\begin{array}{c} t \\ u \\ v \end{array} \right]$.
% _d_ is the ray direction.  Variables _u_ , _v_ are barycentric coordinates 
% and _t/|d|_ is the distance from the intersection point to the ray origin. 
% Ray/triangle intersect if all _t_, _u_, _v_ and _w_=1-u-v are positive.
% # count ray / surface intersections
% # even number means inside and odd mean outside
% # in rare case the ray hits one of the surface faces right on the
%    edge repeat the process with a new ray
%
%% References
% Based on
% * "Fast, minimum storage ray-triangle intersection". Tomas Möller and
%    Ben Trumbore. Journal of Graphics Tools, 2(1):21--28, 1997.
%    http://www.graphics.cornell.edu/pubs/1997/MT97.pdf (Ray/triangle
%    intersection)
% * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/raytri.c
%    (Ray/triangle intersection)
% * Robert Sedgewick "Algorithms" (point in polygon algorithm)
%
%% Licence
% *The function is distributed under BSD License*
format compact; % viewing preference
clear variables; close all;
rng('shuffle');
type('license.txt')

%% Test if random points are inside or outside of the volume
% defined by MATLAB test object "tetmesh"
load tetmesh;
TR = triangulation(tet,X);
[S.faces, S.vertices] = freeBoundary(TR);
n = 2000; % number of points
points = 80*rand(n,3) - repmat([40 40 0], n, 1);
tic
in1 = in_polyhedron(S, points);
fprintf('Number of points inside is %i, outside is %i. Calculation time: %f sec\n', ...
  nnz(in1), nnz(in1==0), toc);

%% Plot results
figure, hold on, view(3)        % Display the result
set(gcf, 'Position', get(gcf, 'Position').*[0 0 1.5 1.5])
patch(S,'FaceColor','g','FaceAlpha',0.2)
plot3(points( in1,1),points( in1,2),points( in1,3),'bo','MarkerFaceColor','b')
plot3(points(~in1,1),points(~in1,2),points(~in1,3),'r.'), axis image
legend({'volume', 'points inside', 'points outside'}, 'Location', 'southoutside')

%% Compare the results to the output of similar inpolyhedron function
% by Sven Holcombe (http://www.mathworks.com/matlabcentral/fileexchange/37856)
% inpolyhedron function is usually faster but requires knowlege about the
% face normals.
if exist('inpolyhedron.m', 'file')
  tic
  in2 = inpolyhedron(S, points);
  fprintf('Number of points inside is %i, outside is %i. Calculation time: %f sec\n', ...
    nnz(in1), nnz(in1==0), toc);
  fprintf('Number of differences is %i\n', sum(in1~=in2));
end

%% Flip 50% of face normals and repeat
msk = rand(size(S.faces,1),1) > 0.5;
S.faces(msk,:) = fliplr(S.faces(msk,:));
in3 = in_polyhedron(S, points);
fprintf('Number of differences for in_polyhedron is %i\n', sum(in1~=in3));
if exist('inpolyhedron.m', 'file')
  in4 = inpolyhedron(S, points);
  fprintf('Number of differences for  inpolyhedron is %i\n', sum(in1~=in4));
end