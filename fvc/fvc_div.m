function RHS = fvc_div(faceTensorField, Sf)
%fvc_div For each face in the mesh it computes the 
% contribution from the divergent discretization [Darwish - 2009].
%   Parameters:
% faceTensorField = Any surfaceTensorField.
% Sf = The surface normal vectors.
%   Return:
% RHS = The Right-hand-side divergent contribution to the final linear
% system. 

%% Initialization
mesh = faceTensorField.mesh;
RHS = zeros(mesh.numberOfElements,3);
fluxField = faceTensorField*Sf;

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
