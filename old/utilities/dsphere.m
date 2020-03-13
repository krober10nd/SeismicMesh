function d=dsphere(p,xc,yc,zc,r)

d=sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2+(p(:,3)-zc).^2)-r;