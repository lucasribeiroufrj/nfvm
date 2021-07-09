function volGradPhi = fvc_GaussGreenGrad(facePhi, Sf)
% fvc_GaussGreenGrad Calculate gradient at cell centroids using 
% compact stencil method (option 2) [Pag. 278 - Moukalled].

%% BEG - Internal faces contribution to the gradient
mesh = facePhi.mesh();

% Internal gradient
iGradPhiData = zeros(3,3,mesh.numberOfElements);

sf = Sf.internalField;
phi = facePhi.internalField;

deltaGrad = phi & sf;

faceOwnerLists = mesh.faceOwnerLists();
for iList = 1:length(faceOwnerLists)
    
    iOwners = faceOwnerLists{iList}.iOwner(:);
    iFaces = faceOwnerLists{iList}.faces(:);
    
    iGradPhiData(:,:,iOwners) = ...
        iGradPhiData(:,:,iOwners) + deltaGrad.range(iFaces);
end

faceNeighLists = mesh.faceNeighLists();
for iList = 1:length(faceNeighLists)
    
    iNeighs = faceNeighLists{iList}.iNeighbour(:);
    iFaces = faceNeighLists{iList}.faces(:);
    
    iGradPhiData(:,:,iNeighs) = ...
        iGradPhiData(:,:,iNeighs) - deltaGrad.range(iFaces);
end
%% END

%% BEG - Boundary faces contribution to the gradient
for iBoundary=1:mesh.numberOfBoundaries
    
    boundary = mesh.boundaries(iBoundary);

    if strcmp(boundary.type,'empty') == true
        continue;
    end

    startFace = boundary.startFace;
    endFace = startFace + boundary.numberOfBFaces - 1;
    owners = [mesh.faces(startFace:endFace).iOwner];
    
    sf = Sf.boundaryField(iBoundary).field;
    phi = facePhi.boundaryField(iBoundary).field;
    
    deltaGrad = phi & sf;
    
    iGradPhiData(:,:,owners) = iGradPhiData(:,:,owners) + deltaGrad.data;
end
%% END

%% Divide by volume to finally obtain the gradient at the volume centroid
iGradPhiField = tensorField(iGradPhiData) / mesh.volume();

%% BEG - Set boundary gradient equal to associated element gradient.
bGradPhiField = bTensorField(mesh.numberOfBoundaries);

for iBoundary=1:mesh.numberOfBoundaries
    
    boundary = mesh.boundaries(iBoundary);
    
    if strcmp(boundary.type,'empty') == true

        pField = patchTensorField(mesh, iBoundary);
    else
        startFace = boundary.startFace;
        endFace = startFace + boundary.numberOfBFaces - 1;

        owners = [mesh.faces(startFace:endFace).iOwner];

        pField = ...
            patchTensorField...
            (...
                mesh,...
                iBoundary,...
                iGradPhiField.subField(owners)...
            );
    end
    
    % Add the calculated patch to Sigma's boundary field.
	bGradPhiField = ...
        bGradPhiField.setPatch...
        (...
            iBoundary,...
            pField...
        );
end
%% END

% "Mount" the final volGradPhi field using the parts: 
% internalfield + boundaryField
volGradPhi = volTensorField(mesh, iGradPhiField, bGradPhiField);

end

