function volGradPhi = fvc_LeastSquareGrad(volPhi)
%fvc_LeastSquareGrad Calculate gradient at cell centroids using 
% Least-Square gradient [Pag. 285 - Moukalled].
%
% Note that the matrix A depends only on the geometry of the mesh, thus if
% the mesh is not deforming, it can be computed only once during the whole
% simualation [muzaferija - 1996].
%
% Ref.:
% [muzaferija - 1996] - Finite-Volume CFD Procedure and Adaptive Error 
% Control Strategy for Grids of Arbitrary Topology
%
% TODO:
% - Check if the simulation does not deform the mesh (see above).

mesh = volPhi.mesh();
numberOfElements = mesh.numberOfElements;
phi = volPhi.internalField().data;
gradPhi = zeros(3,3,numberOfElements);

for iElement=1:numberOfElements
    
    element = mesh.elements(iElement);
    
    for cmpt = 1:3
        
        A = zeros(3, 3);
        B = zeros(3, 1);
        
        for i=1:element.numberOfNeighbours
            
            iNeighbour = element.iNeighbours(i);
            neighbour = mesh.elements(iNeighbour);
            rCF = neighbour.centroid - element.centroid;
            
            w = 1/norm(rCF, 2);
            
            deltaX = rCF(1);
            deltaY = rCF(2);
            deltaZ = rCF(3);
            
            A(1,1) = A(1,1) + w*deltaX*deltaX;
            A(1,2) = A(1,2) + w*deltaX*deltaY;
            A(1,3) = A(1,3) + w*deltaX*deltaZ;
            A(2,1) = A(2,1) + w*deltaY*deltaX;
            A(2,2) = A(2,2) + w*deltaY*deltaY;
            A(2,3) = A(2,3) + w*deltaY*deltaZ;
            A(3,1) = A(3,1) + w*deltaZ*deltaX;
            A(3,2) = A(3,2) + w*deltaZ*deltaY;
            A(3,3) = A(3,3) + w*deltaZ*deltaZ;

            deltaPhi = phi(:,iNeighbour) - phi(:,iElement);
            
            B(1) = B(1) + w*deltaX*deltaPhi(cmpt);
            B(2) = B(2) + w*deltaY*deltaPhi(cmpt);
            B(3) = B(3) + w*deltaZ*deltaPhi(cmpt);
        end

        C = zeros(3,3);
        C(cmpt,:) = A\B;
        gradPhi(:,:,iElement) = gradPhi(:,:,iElement) + C;
    end
end
%%


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
                tensorField(gradPhi(owners))...
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

% "Mount" the final volGradPhi field using the parts: 
% internalfield + boundaryField
volGradPhi = volTensorField(mesh, gradPhi, bGradPhiField);
%% END - Set boundary gradient equal to associated element gradient.

end
