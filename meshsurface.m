function [p,t]=meshsurface(cp,fh,h0,bbox,itmax)
%MESHSURFACE 3-D Surface Mesh Generator using Bezier surfaces.
%   [P,T]=MESHSURFACE(CP,FH,H0,BBOX,FPARAMS)
%    kjr, usp, 2019 

dptol=1e-4; ttol=.1; Fscale=1.2; deltat=.2; deps=sqrt(eps)*h0;

%1a. Create distribution of points in unit plane
pinit{1} = bbox(1,1):h0:bbox(2,1);
pinit{2} = bbox(1,2):h0*sqrt(3)/2:bbox(2,2);
pp=cell(1,2); [pp{:}]  = ndgrid( pinit{:} );
pp{1}(:,2:2:end) = pp{1}(:,2:2:end) + h0/2;
p(:,1) = pp{1}(:); p(:,2) = pp{2}(:);
%1b. Map points to unit plane [0,1] x[0,1]
up=mapToUnitSpace(p) ;
%1c. Project to Bezier surface in R^3 and calculate r0.
[bp(:,1),bp(:,2),bp(:,3)] = BezierSurface(cp,up(:,1),up(:,2));
r0=fh(bp);                                    % Probability to keep point
up=up(rand(size(bp,1),1)<min(r0)^2./r0.^2,:);
up=[0,0; 0,1 ; 1,1; 1,0; up];                 % Add corner points
t = delaunay(up) ; p=[];
[p(:,1),p(:,2),p(:,3)] = BezierSurface(cp,up(:,1),up(:,2));

% Form connectivities (for trisurfupd)
[t2t,t2n]=mkt2t(t);
t2t=int32(t2t-1)'; t2n=int8(t2n-1)';

N=size(p,1);                                         % Number of points N
pold=inf;                                            % For first iteration
it = 1 ; 
disp('Meshing...'); 
while 1
  p0=p;
  % 3. Retriangulation by edge flips 
  if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
    pold=p;                                          % Save current positions
    [t,t2t,t2n]=trisurfupd(int32(t-1)',t2t,t2n,p');  % Update triangles
    t=double(t+1)';
    % 4. Describe each bar by a unique pair of nodes
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
    bars=unique(sort(bars,2),'rows');                % Bars as node pairs
    % 5. Graphical output of the current mesh
    clf,patch('faces',t,'vertices',p,'facecol',[.8,.9,1],'edgecol','k');
    title(['Iteration # ',num2str(it)])
    axis equal;axis off;view(3);cameramenu;drawnow
  end
  
  % 6. Move mesh points based on bar lengths L and forces F
  barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
  L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
  hbars=fh( (p(bars(:,1),:)+p(bars(:,2),:))/2);
  L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
  F=max(L0-L,0);                                     % Bar forces (scalars)
  Fvec=F./L*[1,1,1].*barvec;                         % Bar forces (x,y,z components)
  Ftot=full(sparse(bars(:,[1,1,1,2,2,2]),ones(size(F))*[1,2,3,1,2,3],[Fvec,-Fvec],N,3));
  p=p+deltat*Ftot;                                   % Update node positions in R^3
  
  % 7. Bring all points back to the surface using binary search
  for nn = 1 : length(p)
     [p(nn,1),p(nn,2),p(nn,3),up(nn,1),up(nn,2)]=ApproxProj(cp,p(nn,:));
  end
  
  %8. Project points back in the [0x1] x [0x1] domain
  d=drectangle(up); ix=d>0;                                % Find points outside (d>0)
  dgradx=(drectangle([up(ix,1)+deps,up(ix,2)])-d(ix))/deps; %    Numerical
  dgrady=(drectangle([up(ix,1),up(ix,2)+deps])-d(ix))/deps; %    gradient
  up(ix,:)=up(ix,:)-[d(ix).*dgradx,d(ix).*dgrady];     % Project back to boundary
  [p(:,1),p(:,2),p(:,3)] = BezierSurface(cp,up(:,1),up(:,2));

  it = it + 1 ;
  if mod(it,5)==0
      disp(['Iteration ',num2str(it),' complete']);
  end
  
  % 9. Termination criterion: All nodes move less than dptol (scaled)
  %    or exhausted iterations
  if max(sqrt(sum((p-p0).^2,2))/h0)<dptol, break; end
  if it == itmax, disp('Exhausted iterations'); break; end
end

return; 

% distance function to reproject points that've exited the unit domain 
function d = drectangle(p)
    d =-min(min(min(p(:,2),1-p(:,2)),p(:,1)),1-p(:,1));
