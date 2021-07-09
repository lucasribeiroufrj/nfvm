function [LHS, RHS] = fvm_laplacian(impK, surfaceGradPhi, phi)
%fvm_laplacian For each face in the mesh it computes the 
% coefficients given in (16) and (17) [Darwish - 2009].
%   Parameters:
% surfaceGradPhi = Is the gradient(phi) as a face field.
% Sf = The surface vectors. It does not use mesh.faces.Sf to handle
% anisoptropic diffusion. The mesh.faces.Sf is always untouched, but the
% Sf might have been transformed by the diffusion coefficient Gamma
%   Return:
% LHS = The Left-hand-side laplacian contribution to the final linear
% system. 
% RHS = The Right-hand-side laplacian contribution to the final linear
% system.
%
% ATTENTION:
% - A possible non-conjunctional correction 
% [p. 237, Ferziger - 2001].
%
% Ref.:
% [Ferziger - 2001] - Computational methods for fluid dynamics

%% Initialization.
mesh = surfaceGradPhi.mesh();
LHS = zeros(mesh.numberOfElements,mesh.numberOfElements);
RHS = zeros(mesh.numberOfElements,3);

%% BEG - Calculate the internal face fluxes.
Sf = mesh.Sf().internalField();
nf = mesh.nf().internalField();
eCN = mesh.eCN().internalField();
impk = impK.internalField();

[Ef, Tf] = fvmSfDecomposition(Sf, nf, impk, eCN);

CN = mesh.CN().internalField();
EfNorm = Ef.norm();
gDiff = getData(EfNorm / CN.norm());
gradPhi = surfaceGradPhi.internalField();
FluxVf = getData(gradPhi * Tf);
%% END

%% BEG - Assembles diagonal coefficients and RHS's coeff.
% TODO - Merge both For's (create a method).

faceNeighLists = mesh.faceNeighLists();

for iList = 1:length(faceNeighLists)
    
    iNeighs = faceNeighLists{iList}.iNeighbour(:);
    iFaces = faceNeighLists{iList}.faces(:);
    
    flux = gDiff(iFaces);
    fluxVf = FluxVf(:,iFaces);
    
    % TODO: move it to fvMesh
    % Get linear indices to allow multiple assigns
    idxNeighs = sub2ind(size(LHS), iNeighs, iNeighs);
    
    LHS(idxNeighs) = LHS(idxNeighs) - flux';
    RHS(iNeighs,:) = RHS(iNeighs,:) + fluxVf.';
end

faceOwnerLists = mesh.faceOwnerLists();

for iList = 1:length(faceOwnerLists)
    
    iOwners = faceOwnerLists{iList}.iOwner(:);
    iFaces = faceOwnerLists{iList}.faces(:);
    
    flux = gDiff(iFaces);
    fluxVf = FluxVf(:,iFaces);
    
    % TODO: move it to fvMesh
    % Get linear indices to allow multiple assigns
    idxOwners = sub2ind(size(LHS), iOwners, iOwners);
    
    LHS(idxOwners) = LHS(idxOwners) - flux';
    RHS(iOwners,:) = RHS(iOwners,:) - fluxVf.';
end
%% END

%% BEG - Assembles non-diagonal coefficients.
iOwners = mesh.iOwners();
iNeighs = mesh.iNeighs();
iFaces = 1:mesh.numberOfInteriorFaces;

flux = gDiff(iFaces);

% TODO: move it to fvMesh
% Get linear indices to allow multiple assigns
idxOwnerNeigh = sub2ind(size(LHS), iOwners, iNeighs);
idxNeighOwner = sub2ind(size(LHS), iNeighs, iOwners);

LHS(idxOwnerNeigh) = flux;
LHS(idxNeighOwner) = flux; 
%% END

%% BEG - Add boundary contribution from laplacian operator
for iBoundary=1:mesh.numberOfBoundaries

    boundary = mesh.boundaries(iBoundary);
    
    if strcmp(boundary.type,'empty') == true
        continue;
    end
    
    boundaryField = phi.boundaryField(iBoundary);
    [ fluxLHS, fluxRHS ] = boundaryField.divImpKgradUfluxes();
    owners = boundary.owners;

    for o = 1:length(owners)

        owner = owners(o);

        if ~isempty(fluxLHS)
            
            LHS(owner,owner) = LHS(owner,owner) + fluxLHS(o);
        end
        
        if ~isempty(fluxRHS)
        
            RHS(owner,:) = RHS(owner,:) + fluxRHS(o,:);
        end
    end
end
%% END

end
