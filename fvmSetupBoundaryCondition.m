function valueAtBoundary = fvmSetupBoundaryCondition(...
    theMesh, type, boundary, faceTensorField)
%FVMSETUPBOUNDARYCONDITION Summary of this function goes here
%   Detailed explanation goes here

    function value = calcStress(face, iFace)
        
        sigma = faceTensorField(iFace);
        normal = face.Sf / norm(face.Sf);
        value = sigma * normal;
    end

    function computeBoundaryvalues(callback)
        
        startFace = boundary.startFace;
        endFace = startFace + boundary.numberOfBFaces - 1;

        for iFace=startFace:endFace

            face = theMesh.faces(iFace);
            valueAtBoundary(iFace - startFace + 1) = callback(face, iFace);
        end
    end

valueAtBoundary = vectorField(boundary.numberOfBFaces);

if strcmp(type,'Dirichlet')

    computeBoundaryvalues(@(face, iFace) faceTensorField(iFace));
    
elseif strcmp(type, 'zeroDisplacement')
    
    % Alreary set (zero)
    return;
    
elseif strcmp(type, 'solidTraction')
    
    computeBoundaryvalues(@calcStress);
    
elseif strcmp(type, 'uniformTraction')
    
    computeBoundaryvalues(@(~, ~) faceTensorField )
    
elseif strcmp(type, 'zeroTraction')
    
    % Alreary set (zero)
    return;
    
else
    ME = MException('Error: Unknown boundary condition type: %s',type);
    throw(ME)
end

end

