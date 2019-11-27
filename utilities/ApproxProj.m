function [Xnew,Ynew,Znew,u0,v0] = ApproxProj( cp, qp, w )
% Projects the query point in R^3 onto the Bezier surface defined by the
% control points located in the cell-array cp using a binary search
% algorithm. 
% Kjr, usp, 2019 

thres = 1e-6; % threshold to stop binary search. 

% Initial grid to find starting point 
N=2; 
[u,v]=meshgrid(linspace(0,1,N),linspace(0,1,N)); 

[X,Y,Z]=BezierSurface(cp,u(:),v(:),w);

% compute distance from qp to coarse representation of surface 
for i = 1 : numel(u)
   td(i) = sqrt( (X(i) - qp(1)).^2 + ...
                 (Y(i) - qp(2)).^2 + ...
                 (Z(i) - qp(3)).^2); 
end

minTd = min(td);
id=find(td==minTd,1,'first');
u0=u(id) ; v0=v(id) ; % starting point for binary search

interval = 1/N ; halfInterval = interval/2;

count  = 0 ; 
update = 0 ; 
Tests=[];
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
    % projecting points in the parameter space to the 
    % Bezier surface. 
    [Xr,Yr,Zr]=BezierSurface(cp,Tests(:,1),Tests(:,2),w);
    
    % is there a distance shorter than the current point? 
    for i = 1 : 8
        tdd(i) = sqrt( (Xr(i) - qp(1)).^2 + ...
            (Yr(i) - qp(2)).^2 + ...
            (Zr(i) - qp(3)).^2);
    end
    
    minTdd = min(tdd) ;
    idd = find(tdd==minTdd,1,'first');

    % all the distances are larger. Lets try again but halve the search space. 
    if(update || minTdd >= minTd)
        interval = interval/2;
        halfInterval = halfInterval/2;
    end
    
    % all the distances are smaller, lets start from this position
    if(minTdd < minTd)
        update=1 ; 
        u0 = Tests(idd,1);
        v0 = Tests(idd,2);
    end
    
    count = count + 1;
    
    if interval < thres
        break
    end
end
[Xnew,Ynew,Znew]=BezierSurface(cp,u0,v0,w);


