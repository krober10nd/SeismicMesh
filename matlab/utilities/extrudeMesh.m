function [pAbove,pBelow]=extrudeMesh(p,t)
%% Create nodes above and below mesh facets at a distance of the average 
%  edgelength. 
tria = triangulation(t,p);
fNvec=tria.faceNormal;
for i = 1 : length(t)
    t1=t(i,1); t2=t(i,2); t3=t(i,3);
    for j = 1 : 3
        pmid(i,:) = (p(t1,:) + p(t2,:) + p(t3,:))./3;
    end
end
for i =1 : length(t)
    bars=[t(i,[1,2]);t(i,[1,3]);t(i,[2,3])];
    barvec=p(bars(:,1),:)-p(bars(:,2),:);
    L=sqrt(sum(barvec.^2,2));
    meanBarvec(i,1) = (L(1)+L(2)+L(3))/3;
    pAbove(i,:) = pmid(i,:) + fNvec(i,:)*meanBarvec(i);
    pBelow(i,:) = pmid(i,:) - fNvec(i,:)*meanBarvec(i);
end