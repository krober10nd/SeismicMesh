function [Xout,Yout,Zout] = BezierSurface(cp,U,V)
%
% Fast calculation of the point cloud of a Bezier surface
%
% [Xout,Yout,Zout] = BezierSurface(cp,U,V)
%
% Input parameters
%     cp ({n x m}) is a cell containing the control points of the
%         Bezier surface.
%     U ([1 x nU]) is a vector containing the values in the 1-direction
%         (the first dimension of the matrix of points (x,y) at which the
%         Bezier surface is plotted). max(U) must be <=1.
%     V ([1 x nV]) is a vector containing the values in the 2-direction
%         (the second dimension of the matrix of points (x,y) at which the
%         Bezier surface is plotted). max(V) must be <=1.
%
% Output parameters
%     Xout ([nU x nV]) is the X-coordinates of the points of the
%         Bezier surface.
%     Yout ([nU x nV]) is the Y-coordinates of the points of the
%         Bezier surface.
%     Zout ([nU x nV]) is the Z-coordinates of the points of the
%         Bezier surface.
%
% kjr, usp, 2019

%% Initial checks
[n,m]=size(cp);
n=n-1; % 1-direction
m=m-1; % 2-direction

% convert plotting vectors into matrices
p2=numel(U);

X=zeros(p2,(n+1)*(m+1)); 
Y=zeros(p2,(n+1)*(m+1)); 
Z=zeros(p2,(n+1)*(m+1)); 

k=1;
for i=0:n
    for j=0:m
        niF=factorial(n)/(factorial(i)*factorial(n-i));
        Bin=niF*V.^i.*(1-V).^(n-i);
        mjF=factorial(m)/(factorial(j)*factorial(m-j));
        Bjm=mjF*U.^j.*(1-U).^(m-j);
        X(:,k)=Bin.*Bjm.*cp{i+1,j+1}(1);
        Y(:,k)=Bin.*Bjm.*cp{i+1,j+1}(2);
        Z(:,k)=Bin.*Bjm.*cp{i+1,j+1}(3);
        k=k+1;
    end
end
Xout=sum(X,2); 
Yout=sum(Y,2); 
Zout=sum(Z,2); 

end

