function [pmap,mn,mx]=mapToUnitSpace(p)
% map coordinates to unit space
mn=min(p);

p = p - repmat(mn,[size(p,1),1]);

mx=max(p);

pmap = p./repmat(mx,[size(p,1),1]);

end