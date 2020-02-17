function inside = in_polyhedron(varargin)
%% Point within the volume test
%IN_POLYHEDRON  Tests if points are on or inside a 3D triangulated surface 
% (faces/vertices) or volume (tetrahedrals/vertices). There is NO
% assumption about orientation of the face normals.
%
%   IN = IN_POLYHEDRON(X,POINTS) tests if the query points (POINTS) are inside
%   the surface/polyhedron defined by X. X can be a structure with fields 
%   'vertices' and 'faces' or an object of MATLAB triangulation class. In
%   case of triangulation class object we will only use the outside
%   boundary. POINTS is an N-by-3 set of XYZ coordinates. IN is an N-by-1 
%   logical vector which will be TRUE for each query point inside the surface.
%
%   IN_POLYHEDRON(FACES,VERTICES,POINTS) takes faces/vertices separately, 
%   rather than in an FV structure.
%
% Algorithm:
% For each point do:
% 1) shoot a random ray out of the query point in a random direction
% 2) for each face solve 
%        |t|
%    M * |u| = (o-v0)
%        |v|
%  for [t; u; v] where M = [-d, v1-v0, v2-v0]. "d" is the ray direction. 
%  u, v, w (=1-u-v) are barycentric coordinates and t is the distance from 
%  the ray origin in units of |d|.
%  Ray/triangle intersect if all t, u, v and w are positive.
% 3) count ray / surface intersections
% 4) even number means inside and odd mean outside
% 5) in rare case the ray hits one of the surface faces right on the
%    edge repeat the process with a new ray
%
% Based on:
%  * Ray/triangle intersection
%    * "Fast, minimum storage ray-triangle intersection". Tomas Möller and
%      Ben Trumbore. Journal of Graphics Tools, 2(1):21--28, 1997.
%      http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
%    * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/raytri.c
%  * point in polygon algorithm: Robert Sedgewick "Algorithms"
%
% Author:
%    Jarek Tuszynski (jaroslaw.w.tuszynski@leidos.com)
%
% License: BSD license (http://en.wikipedia.org/wiki/BSD_licenses)

%% Get POINTS, FACES and VERTICES inputs
if  isa(varargin{1}, 'triangulation') % in_polyhedron(triangulation_class, POINTS)
  [faces, vertices] = freeBoundary(varargin{1});
  points = varargin{2};
elseif isstruct(varargin{1})   % in_polyhedron(FVstruct, POINTS)
  ok = isfield(varargin{1}, 'vertices') && isfield(varargin{1}, 'faces');
  assert(ok, 'Structure FV must have "faces" and "vertices" fields' );
  faces    = varargin{1}.faces;
  vertices = varargin{1}.vertices;
  points   = varargin{2};
else                       % in_polyhedron(FACES, VERTICES, POINTS)
  faces    = varargin{1};
  vertices = varargin{2};
  points   = varargin{3};
end
clear varargin

%% Transpose inputs if needed
if (size(points  ,1)==3 && size(points  ,2)~=3), points   = points';   end
if (size(vertices,1)==3 && size(vertices,2)~=3), vertices = vertices'; end
if (size(faces   ,1)==3 && size(faces   ,2)~=3), faces    = faces';    end
assert(size(points  ,2)==3, '"Points" array must be in Nx3 format');
assert(size(vertices,2)==3, '"Vertices" array must be in Nx3 format');
assert(size(faces   ,2)==3, '"Faces" array must be in Nx3 format');

%% Convert the polyhedron into array of faces, stored as a point and 2 vectors
vert0  = vertices(faces(:,1),:);
edge1  = vertices(faces(:,2),:)-vert0;  % find vectors for two edges sharing vert0
edge2  = vertices(faces(:,3),:)-vert0;
N      = size(vert0,1);
clear vertices faces                    % those are no longer needed

%% In case of 3D data use the following algorithm:
% 1) shoot a random ray out of the point in a random direction
% 2) solve for each face
%        |t|
%    M * |u| = (o-v0)
%        |v|
%  for [t; u; v] where M = [-d, v1-v0, v2-v0]. u,v are barycentric coordinates
%  and t - the distance from the ray origin in |d| units
%  ray/triangle intersect if t>=0, u>=0, v>=0 and u+v<=1
% 3) count ray / surface intersections
% 4) even number means inside and odd mean outside
% 5) in rare case the ray hits one of the surface faces right on the
%    edge repeat the process with a new ray

% distance (in barycentric units (0-1)) from the edge where we might get 
% roundoff errors.
eps    = 1e-10; 
inside = nan+zeros(size(points,1),1); % nan indicates that there was no succesful test yet
while any(isnan(inside))
  dir  = repmat(rand(1,3)-0.5,N,1);   % pick random direction for the ray
  pvec = cross_prod(dir, edge2);
  det  = sum(edge1.*pvec,2); % determinant of the matrix M = dot(edge1,pvec)
  angleOK = (abs(det)>eps);  % if determinant is near zero then ray lies in the plane of the triangle
  
  %% For each point which we did not calculated yet...
  for iPoint = 1:size(points,1)
    if ~isnan(inside(iPoint))
      continue
    end
    tvec = -bsxfun(@minus, vert0, points(iPoint,:)); % vector from vert0 to ray origin
    u    = sum(tvec.*pvec,2)./det;    % 1st barycentric coordinate
    
    % limit some calculations only to line/triangle pairs where it makes
    % a difference. It is tempting to try to push this concept of
    % limiting the number of calculations to only the necessary to "u"
    % and "t" but that produces slower code
    ok = (angleOK & u>-eps & u<=1.0+eps); % mask out the faces that AFAIK do not intersect
    u = u(ok,:); % trim so it is the size of v and t
    % if all line/plane intersections are outside the triangle than no intersections
    if ~any(ok)
      if ~any(u>eps & u<=1.0-eps) % if far away from the edges...
        inside (iPoint) = false;
      end
      continue
    end
    qvec = cross_prod(tvec(ok,:), edge1(ok,:)); % prepare to test V parameter
    v = sum(dir  (ok,:).*qvec,2)./det(ok,:); % 2nd barycentric coordinate
    t = sum(edge2(ok,:).*qvec,2)./det(ok,:); % distance from point to triangle
    % test if line/plane intersection is within the triangle
    bary = [u, v, 1-u-v, t];      % barymatric coordinates
    intersect = all(bary>-eps,2); % intersection on the correct side of the origin
    baryi = bary(intersect,:);    % barymatric coordinates for the intersections only
    if all( min(abs(baryi),[], 2)>eps ) % ray did not hit the edge of the triangle
      nIntersect = sum(intersect);    % number of intersections
      inside(iPoint) = mod(nIntersect,2)>0; % inside if odd number of intersections
    end
    %if any(abs(max(bary,[], 2)-1)<eps & abs(t)<eps) % test for point being one of the grid points
    if any(max(bary,[], 2)<1+eps & abs(t)<eps) % test for point belonging to the surface
      inside(iPoint) = true;
    end
  end
end
inside = (inside~=0); % convert to boolean

%% ========================================================================
function c=cross_prod(a,b)
% strip down version of MATLAB cross function
c = [a(:,2).*b(:,3)-a(:,3).*b(:,2), a(:,3).*b(:,1)-a(:,1).*b(:,3), a(:,1).*b(:,2)-a(:,2).*b(:,1)];
