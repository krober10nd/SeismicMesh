function [Xnew,Ynew,Znew,u0,v0] = ApproxProj( cp, qp, w, ustart,vstart )
% Projects the query point in R^3 onto the Bezier surface defined by the
% control points located in the cell-array cp using a binary search
% algorithm. 
% Kjr, usp, 2019 

thres = 1e-6; % threshold to stop binary search. 

u0 = ustart; v0 = vstart;

[X,Y,Z]=BezierSurface(cp,u0,v0,w);

td=(X - qp(1)).^2+(Y-qp(2)).^2 +(Z-qp(3)).^2; 
minTd = min(td);

interval = 1e-3; halfInterval = interval/2;
count  = 0 ; 
update = 0 ; 
Tests=[];
%disp(['min dis at begin is ',num2str(sqrt(minTd))]); 
while 1
    
    % grid around the point in question 
    Tests(1,:) = [u0+halfInterval, v0];
    Tests(2,:) = [u0+halfInterval, v0+halfInterval];
    Tests(3,:) = [u0, v0+halfInterval];
    Tests(4,:) = [u0-halfInterval, v0+halfInterval];
    Tests(5,:) = [u0-halfInterval, v0];
    Tests(6,:) = [u0-halfInterval, v0-halfInterval];
    Tests(7,:) = [u0, v0-halfInterval];
    Tests(8,:) = [u0+halfInterval,v0-halfInterval];
    
    % 75% of all the time is spent in this line
    [Xr,Yr,Zr]=BezierSurface(cp,Tests(:,1),Tests(:,2),w);
    
    % is there a distance shorter than the current point?
    tdd=(Xr-qp(1)).^2 +(Yr-qp(2)).^2 +(Zr-qp(3)).^2;
    minTdd = min(tdd) ;
    idd = find(tdd==minTdd,1,'first');
    
    % all the distances are larger. Lets try again but halve the search space.
    if(update || minTdd >= minTd)
        interval = interval/2;
        halfInterval = halfInterval/2;
    end
    
    % all the distances are smaller, lets start from this position
    if(minTdd < minTd)
        minTd=minTdd;
        update=1 ; 
        u0 = Tests(idd,1);
        v0 = Tests(idd,2);
    end
    
    count = count + 1;
    
    %saveMinTd(count) = minTd;
    
    if interval < thres
        %disp(['min dis at begin is ',num2str(sqrt(minTd))]);
        break
    end
end
[Xnew,Ynew,Znew]=BezierSurface(cp,u0,v0,w);


