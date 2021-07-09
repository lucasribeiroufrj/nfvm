function plotDisplacement(mesh, U, nDivisions, U_analyticalHandle)
%plotDisplacement Summary of this function goes here
%   Detailed explanation goes here

analytical = false;
if nargin == 4
    analytical = true;
end

if analytical
    U = zeros(2,mesh.numberOfElements);
end

X = zeros(1,mesh.numberOfElements);
Y = zeros(1,mesh.numberOfElements);

for iElement=1:mesh.numberOfElements
    x = mesh.elements(iElement).centroid(1);
    y = mesh.elements(iElement).centroid(2);

    if analytical

        [u, v] = U_analyticalHandle(x,y);
        U(1:2,iElement) = [u; v];
    end

    X(iElement) = x;
    Y(iElement) = y;

    u = U(1,iElement);
    v = U(2,iElement);
end

detlaX = (max(X) - min(X))*0.01;
detlaY = (max(Y) - min(Y))*0.01;
intervalX = [ min(X)-detlaX max(X)+detlaX ];
intervalY = [ min(Y)-detlaY max(Y)+detlaY ];
dx = (intervalX(2) - intervalX(1)) / nDivisions;
dy = (intervalY(2) - intervalY(1)) / nDivisions;
[xq,yq] = meshgrid(...
    intervalX(1):dx:intervalX(2),...
    intervalY(1):dy:intervalY(2));
vq = griddata(X,Y,U(1,:),xq,yq);
uq = griddata(X,Y,U(2,:),xq,yq);

%if isempty(filter) == false && false
%    idxToRemove = filter(X,Y);
%    Xg(idxToRemove) = NaN;
%    Yg(idxToRemove) = NaN;    
%end

axis equal;
quiver(xq,yq,vq,uq,1);
detlaX = 5*detlaX;
detlaY = 5*detlaY;
intervalX = [ min(X)-detlaX max(X)+detlaX ];
intervalY = [ min(Y)-detlaY max(Y)+detlaY ];
axis( [ intervalX intervalY ] );
axis square;

end
