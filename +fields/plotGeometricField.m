function quiverObject = plotGeometricField(field, params)
% plotGeometricField Plot a geometricField.
%   quiverObject = plotGeometricField(field, params) plots the 
%   geometricField field and returns the quiver object returned by quiver3
%   call. Use quiverObject.Visible = false; to make the plot invisible.

if nargin == 1
    params = [];
end

nElements = field.mesh().numberOfElements();
A = zeros(6, nElements);

shouldCreateFig = tryGetOrDefault(params, 'createFigure', false);
shouldHoldOn = tryGetOrDefault(params, 'holdOn', true);
autoScale = tryGetOrDefault(params, 'AutoScale', false);
alphaValue = tryGetOrDefault(params, 'alpha', 0.05);

if shouldCreateFig
    figure;
end

if shouldHoldOn
    hold on;
end

centroids = field.mesh().r_C().internalField();
vectors = field.internalField();
A(1:3,:) = centroids.data;
A(4:6,:) = vectors.data;

quiverObject = quiver3...
    (...
        A(1,:), A(2,:), A(3,:), A(4,:), A(5,:), A(6,:),...
        'LineWidth', 1.2,...
        'AutoScale', autoScale ...
    );
alpha(alphaValue);

if shouldHoldOn
    hold off
end

end

