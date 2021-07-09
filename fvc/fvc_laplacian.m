function RHS = fvc_laplacian(impK, faceGradPhi, Sf)
%fvc_laplacian For each face in the mesh it computes the 
% coefficients given in (16) and (17) [Darwish - 2009]. It handles Handles
% anisotropic diffusion [Darwish - 2009].
%   Parameters:
% faceGradPhi = Is the gradient(phi) as a face field.
% Sf = The surface normal vectors.
%   Return:
% RHS = The Right-hand-side laplacian contribution to the final linear
% system. 

%% Initialization
mesh = faceGradPhi.mesh();
RHS = zeros(mesh.numberOfElements,3);
fluxField = faceGradPhi * (impK.' * Sf);

%% BEG - Add contribution from internal faces
faceNeighLists = mesh.faceNeighLists();

for iList = 1:length(faceNeighLists)
    
    iNeighs = faceNeighLists{iList}.iNeighbour(:);
    iFaces = faceNeighLists{iList}.faces(:);
    
    flux = fluxField.internalField(iFaces);
    
    RHS(iNeighs,:) = RHS(iNeighs,:) + flux';
end

faceOwnerLists = mesh.faceOwnerLists();

for iList = 1:length(faceOwnerLists)
    
    iOwners = faceOwnerLists{iList}.iOwner(:);
    iFaces = faceOwnerLists{iList}.faces(:);
    
    flux = fluxField.internalField(iFaces);
    
    RHS(iOwners,:) = RHS(iOwners,:) - flux';
end
%% END

%% BEG - Add contribution from boundary faces
for iBoundary=1:mesh.numberOfBoundaries

    boundary = mesh.boundaries(iBoundary);
    
    if strcmp(boundary.type,'empty') == true
        continue;
    end
 
    iOwners = boundary.owners;
    flux = fluxField.boundaryField(iBoundary).field().data';
    
    for o = 1:length(iOwners)
        
        iOwner = iOwners(o);
        
        RHS(iOwner,:) = RHS(iOwner,:) - flux(o,:);
    end
end
%% END


end
