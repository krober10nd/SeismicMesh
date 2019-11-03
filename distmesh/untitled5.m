bbox = [-2 1
    -2 1
    -2 1];

drectangle3D=@(p,x1,x2,y1,y2,z1,z2)-min(min(min(min(min(-z1+p(:,3),z2-p(:,3)),-y1+p(:,2)),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

fd = @(p) drectangle3D(p,...
    bbox(1,1),bbox(1,2),...
    bbox(2,1),bbox(2,2),...
    bbox(3,1),bbox(3,2));

fh = @(p) ones(size(p,1),1)*.10; 

[p,t]=distmeshnd(fd,fh,0.10,bbox',[]);

figure; simpplot(p,t)