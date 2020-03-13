function [poly,poly_idx]=IdentifyBoundaries(P,T,DIR)

[bnde,bpts]=extdom_edges2(T,P);

% use this to figure out the vstart and vend
figure, plot(bpts(:,1),bpts(:,2),'k.');
%hold on; fastscatter(obj.p(:,1),obj.p(:,2),obj.b) ;
caxis([-10 10]) ; axis equal ;
title('use data cursor to identify vstart and vend');
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn2,bpts})

% Use data cursor to get the values of the boundary nodes.
vstart=input('Enter the value for vstart : ');
vend=input('Enter the value for vend : ');
bndidx=unique(bnde(:));

vstart= bndidx(vstart);
vend  = bndidx(vend);

[poly,poly_idx] = extract_boundary(vstart,vend,bnde,P,DIR);

