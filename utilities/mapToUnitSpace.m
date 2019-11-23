function [pmap,mn,mx]=mapToUnitSpace(p)
% map coordinates to unit space
mn=min(p);

p = p - mn;

mx=max(p);

pmap = p./mx;

end