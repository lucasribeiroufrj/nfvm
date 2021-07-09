function pointFaces = calculatePointFaces(mesh, iBoundary)
%calculatePointFaces Point face list.
%   calculatePointFaces(mesh, boundary) return an array of structures each
%   having a member called iFaces, which is a list of neighboring faces.

boundary = mesh.boundaries(iBoundary);
iNodes = boundary.iNodes;
nNodes = length(iNodes);
pointFaces(nNodes) = struct('iFaces',[]);

for n = 1:nNodes

    iAllFaces = mesh.nodes(iNodes(n)).iFaces;

    % Only faces from the boundary
    pointFaces(n).iFaces = ...
        intersect(iAllFaces, boundary.iFaces);
end

end
