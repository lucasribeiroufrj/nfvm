function pointNormals = calculatePointNormals(mesh, iBoundary)
%calculatePointNormals Point normal list.
%   calculatePointNormals(mesh, boundary) return an array 3 x nNodes,
%   nNodes is the number boundary nodes. The column j is the point normal.
%   This is calculated using an average.

boundary = mesh.boundaries(iBoundary);
iNodes = boundary.iNodes;
nNodes = length(iNodes);
pointNormals = zeros(3,nNodes);

if ~isfield(boundary, 'pointFaces')
    boundary.pointFaces = calculatePointFaces(mesh, iBoundary);
end

for n = 1:nNodes
    
    iFaces = boundary.pointFaces(n).iFaces;

    for iFace = iFaces

        nf = mesh.faces(iFace).nf;
        pointNormals(:,n) = pointNormals(:,n) + nf;
    end
end

for i = 1:length(pointNormals)

    nf = pointNormals(:,i);
    pointNormals(:,i) = nf / norm(nf);
end

end
