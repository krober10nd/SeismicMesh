function txt = myupdatefcn2(~,event_obj,myarray)
pos = get(event_obj,'Position');
ind = find(abs(myarray(:,1)-pos(1))<eps & abs(myarray(:,2)-pos(2))<eps);
txt = {['X: ',num2str(pos(1))],...
    ['Y: ',num2str(pos(2))],...
    ['Index: ',num2str(ind')]};
end