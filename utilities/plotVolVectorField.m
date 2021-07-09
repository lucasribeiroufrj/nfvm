function plotVolVectorField(field, params)
%plotVolVectorField Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    params = [];
end

nElements = field.mesh().numberOfElements();
x = zeros(1, nElements);
y = zeros(1, nElements);
z = zeros(1, nElements);
u = zeros(1, nElements);
v = zeros(1, nElements);
w = zeros(1, nElements);

centroids = field.mesh().r_C().internalField();
values = field.internalField();

for iElem = 1:nElements
    
    vec = values(iElem);
    c = centroids(iElem);
    
    x(iElem) = c(1);
    y(iElem) = c(2);
    z(iElem) = c(3);
    u(iElem) = vec(1);
    v(iElem) = vec(2);
    w(iElem) = vec(3);
end

shouldCreateFig = tryGetOrDefault(params, 'createFigure', false);
shouldHoldOn = tryGetOrDefault(params, 'holdOn', true);
autoScale = tryGetOrDefault(params, 'AutoScale', false);

if shouldCreateFig
    figure;
end

if shouldHoldOn
    hold on;
end

alpha 0.05;
ax = quiver3(x,y,z,u,v,w,'LineWidth',1.2, 'AutoScale', autoScale);

if shouldHoldOn
    hold off
end

end

