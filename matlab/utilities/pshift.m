function p=pshift(p,x0,y0,z0)

p(:,1)=p(:,1)+x0*p(:,1);
p(:,2)=p(:,2)+y0;
p(:,3)=p(:,3)+z0;

