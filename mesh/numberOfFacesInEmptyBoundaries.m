function nFaces = numberOfFacesInEmptyBoundaries( mesh )
% numberOfFacesInEmptyBoundaries(mesh) Return the total number of faces 
%   beloging to all empty boundaries.

nFaces = 0;

for boundary = mesh.boundaries

    if strcmp(boundary.type,'empty')

        nFaces = nFaces + boundary.numberOfBFaces;
    end
end

end

