function plotStress...
    (...
        volStressField,...
        nDivisions,...
        stressHandle,...
        stressType,...
        contourLines...
    )
%plotStress Summary of this function goes here
%   Detailed explanation goes here

% TODO - This function can be merged into the function 
% plotDisplacement.

mesh = volStressField.mesh;

analytical = false;
if ~isempty(stressHandle)
    analytical = true;
end

if analytical
    W = zeros(1,mesh.numberOfElements);
    %phi = zeros(2,mesh.numberOfElements);
end

X = zeros(1,mesh.numberOfElements);
Y = zeros(1,mesh.numberOfElements);

data = volStressField.internalField().data;

for iElement=1:mesh.numberOfElements
    
    x = mesh.elements(iElement).centroid(1);
    y = mesh.elements(iElement).centroid(2);
    a = stressType(1);
    b = stressType(2);

    if analytical
        
        stress = stressHandle(x,y);
        W(1,iElement) = stress(a,b);
    else
        stress = data(:,:,iElement); %volStressField.internalField(iElement);
        W(1,iElement) = stress(a,b);
    end

    X(iElement) = x;
    Y(iElement) = y;
end

N = nDivisions;
detlaX = (max(X) - min(X))*0.01;
detlaY = (max(Y) - min(Y))*0.01;
intervalX = [ min(X)-detlaX max(X)+detlaX ];
intervalY = [ min(Y)-detlaY max(Y)+detlaY ];
dx = (intervalX(2) - intervalX(1)) / N;
dy = (intervalY(2) - intervalY(1)) / N;
[xq,yq] = meshgrid(...
intervalX(1):dx:intervalX(2),...
intervalY(1):dy:intervalY(2));
vq = griddata(X,Y,W,xq,yq);
%vq = griddata(X,Y,phi(1,:),xq,yq);
%uq = griddata(X,Y,phi(2,:),xq,yq);

[C, h ]= contour(xq, yq, vq, contourLines);%, 'ShowText','on');
set(h,'LineColor','black')
axis equal;
axis square;

end

