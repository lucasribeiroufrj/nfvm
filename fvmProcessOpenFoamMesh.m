function mesh = fvmProcessOpenFoamMesh(mesh,Sf_decomposition)

% TODO - This is terrible!
over_relaxed = 'over_relaxed';

if nargin < 2
    Sf_decomposition = over_relaxed;
end

%% BEG - setup Node connectivities (part 1)
%
for iFace=1:mesh.numberOfFaces
    
    iNodes = mesh.faces(iFace).iNodes;
    for iNode=iNodes
       mesh.nodes(iNode).iFaces = [mesh.nodes(iNode).iFaces iFace];
    end
    
end
%% END - setup Node connectivities (part 1)


%% BEG - setup Node connectivities (part 2)
%
for iElement=1:mesh.numberOfElements
    
    iFaces = mesh.elements(iElement).iFaces;
    
    for iFace=iFaces
        
       iNodes = mesh.faces(iFace).iNodes;
       
       for iNode=iNodes
           
           if(sum(mesh.elements(iElement).iNodes == iNode)==0)
               
              mesh.elements(iElement).iNodes = ...
                  [mesh.elements(iElement).iNodes iNode];
              
              mesh.nodes(iNode).iElements = ...
                  [mesh.nodes(iNode).iElements iElement];
           end
       end
    end
end
%% END - setup Node connectivities (part 2)


%% BEG - Process basic Face Geometry
%
numberOfFaces = mesh.numberOfFaces;
for iFace=1:numberOfFaces
    
    mesh.faces(iFace).edges = [];
    mesh.faces(iFace).DeltaVol=0;
    iNodes = mesh.faces(iFace).iNodes;
    numberOfiNodes = length(iNodes);
    %
    % Compute a rough centre of the face
    %
    centre = [0 0 0]';
    for iNode=iNodes
        centre = centre + mesh.nodes(iNode).centroid;
    end
    centre = centre/numberOfiNodes;

    centroid = [0 0 0]';
    Sf = [0 0 0]';
    area = 0;
    %
    % using the centre compute the area and centoird of vitual triangles
    % based on the centre and the face nodes
    %
    for iTriangle=1:numberOfiNodes
        point1 = centre;
        iNode2 = iNodes(iTriangle);
        point2 = mesh.nodes(iNode2).centroid;
        if(iTriangle<numberOfiNodes)
            iNode3 = iNodes(iTriangle+1);
            point3 = mesh.nodes(iNode3).centroid;
        else
            iNode3 = iNodes(1);
            point3 = mesh.nodes(iNode3).centroid;
        end
        local_centroid = (point1+point2+point3)/3; 
        local_Sf  = 0.5*cross(point2-point1,point3-point1);
        local_area = fvmMagnitude(local_Sf);
 
        centroid = centroid + local_area*local_centroid;
        Sf = Sf + local_Sf;
        area = area + local_area;
        
        edge.iNodes = [ iNode2 iNode3 ];
        edge.vector = point3 - point2;
        edge.norm = norm(edge.vector);
        edge.biNormal = cross(edge.vector/edge.norm, Sf/norm(Sf));
        mesh.faces(iFace).edges = [ mesh.faces(iFace).edges  edge ];
    end
    centroid = centroid/area;
 
    %
    mesh.faces(iFace).centroid = centroid;
    mesh.faces(iFace).Sf = Sf;
    mesh.faces(iFace).nf = Sf / norm(Sf);
    mesh.faces(iFace).area = area;
end
%% END - Process basic Face Geometry


%% BEG - Compute volume and centroid of each element
%
numberOfElements = mesh.numberOfElements;
for iElement=1:numberOfElements
    
    iFaces = mesh.elements(iElement).iFaces;
    %
    % Compute a rough centre of the element
    %
    centre = [0 0 0]';
    for iFace=1:length(iFaces)
        centre = centre + mesh.faces(iFace).centroid;
    end
    centroid = [0 0 0]';
    Sf = [0 0 0]';
    centre = centre/length(iFaces);
    % using the centre compute the area and centoird of vitual triangles
    % based on the centre and the face nodes
    %
    localVolumeCentroidSum = [0 0 0]';
    localVolumeSum = 0;
    for iFace=1:length(iFaces)
        localFace = mesh.faces(iFaces(iFace));
        localFaceSign = mesh.elements(iElement).faceSign(iFace);
        Sf = localFace.Sf*localFaceSign;
        Cf = localFace.centroid - centre;
        % calculate face-pyramid volume
        localVolume = Sf'*Cf/3;
        % Calculate face-pyramid centre
        localCentroid = 0.75*localFace.centroid + 0.25*centre;

        %Accumulate volume-weighted face-pyramid centre
        localVolumeCentroidSum = localVolumeCentroidSum + localCentroid*localVolume;

         % Accumulate face-pyramid volume
        localVolumeSum = localVolumeSum + localVolume;
    end
    centroid = localVolumeCentroidSum/localVolumeSum; 
    volume = localVolumeSum;
    %
    mesh.elements(iElement).volume = volume;
    mesh.elements(iElement).OldVolume = volume;
    mesh.elements(iElement).centroid = centroid;
end
%% END - Compute volume and centroid of each element


%% BEG - Process secondary Face Geometry (internal faces)
%
numberOfInteriorFaces = mesh.numberOfInteriorFaces;
for iFace=1:numberOfInteriorFaces
    %
    theFace = mesh.faces(iFace);
    nf = theFace.Sf/theFace.area;
    %
    element1 = mesh.elements(theFace.iOwner);
    element2 = mesh.elements(theFace.iNeighbour);
    %
    CN = element2.centroid - element1.centroid;
    mesh.faces(iFace).CN = CN;
    
    eCN = CN/fvmMagnitude(CN);
    mesh.faces(iFace).eCN = eCN;
    
    E = theFace.area*eCN;
    if strcmp(Sf_decomposition, over_relaxed) == true
        cosTheta = eCN'*nf;
        E = 1/cosTheta*E;
    end
    
    mesh.faces(iFace).gDiff = fvmMagnitude(E)/fvmMagnitude(CN);
    mesh.faces(iFace).T = theFace.Sf - E;
    mesh.faces(iFace).E = E;

    %
    Cf = theFace.centroid - element1.centroid;
    fF = element2.centroid - theFace.centroid;
    dFC = element2.centroid - element1.centroid;
    mesh.faces(iFace).gf = (Cf'*nf)/(Cf'*nf + fF'*nf);
    mesh.faces(iFace).gc = norm(fF) / norm(dFC);
    mesh.faces(iFace).walldist = 0;
    %
end
%% END - Process secondary Face Geometry


%% BEG - Process secondary Face Geometry (boundary faces)
%
for iFace=numberOfInteriorFaces+1:numberOfFaces
     %
    theBFace = mesh.faces(iFace);
    nf = theBFace.Sf/theBFace.area;
    %
    element1 = mesh.elements(theBFace.iOwner);
    %
    CN = theBFace.centroid - element1.centroid;
    mesh.faces(iFace).CN = CN;
    mesh.faces(iFace).gDiff = theBFace.area* theBFace.area/dot(CN,theBFace.Sf);
    
    eCN = CN/fvmMagnitude(CN);
    mesh.faces(iFace).eCN = eCN;
    
    E = theBFace.area*eCN;
    
    if strcmp(Sf_decomposition, over_relaxed) == true
        cosTheta = eCN'*nf;
        E = 1/cosTheta*E;
    end
    
    mesh.faces(iFace).T = theBFace.Sf - E;
    mesh.faces(iFace).E = E;
    %
    mesh.faces(iFace).gf = 1;
    mesh.faces(iFace).gc = 1;
    mesh.faces(iFace).walldist = (CN'*theBFace.Sf)/norm(theBFace.Sf);
  
end
%% END - Process secondary Face Geometry (boundary faces)


%% BEG - Loop over all the nodes and create a flag to sort boundary nodes
% from interior nodes
%
for iFace = 1:mesh.numberOfInteriorFaces
    
    mesh.faces(iFace).patchIndex = 0;
    iNodes = mesh.faces(iFace).iNodes;
    
    for iNode = iNodes
        
        mesh.nodes(iNode).isInterior=true;
    end 
end

for iBoundary = 1:mesh.numberOfBoundaries 
    
    startFace = mesh.boundaries(iBoundary).startFace;
    endFace = mesh.boundaries(iBoundary).endFace;
    
    for iFace = startFace:endFace

        mesh.faces(iFace).patchIndex = iBoundary;
        iNodes = mesh.faces(iFace).iNodes;

        for iNode=iNodes

            mesh.nodes(iNode).isInterior=false;
        end 
    end
end
%% END - Loop over all the nodes and create a flag to sort boundary nodes
% from interior nodes


%% BEG - Capture boundary elements(faces)
%
% Each node has to know which boundary element* it is connected to. So here
% we create two kind of indices for the same face. One (iBElements) helps
% to mount the global linear system matrix A and the other (iBFaces) is the
% common face indexing.
%
% * Boundary face which does not belong to an empty patch.
%
% Forcing initialization of all nodes.
mesh.nodes(1).iBElements = []; % Indexing relative to cells
mesh.nodes(1).iBFaces = [];    % Indexing relative to faces

for boundary = mesh.boundaries

    if strcmp(boundary.type,'empty')
        continue;
    end
    
    startFace = boundary.startFace;
    endFace = boundary.endFace;
    
    for face = mesh.faces(startFace:endFace)
        
        iFace = face.index;
        iBElement = iFace - mesh.numberOfInteriorFaces + ...
                mesh.numberOfElements;
        mesh.faces(iFace).iBIndex = iBElement;
        
        for n = 1:length(face.iNodes)
            
            iNode = face.iNodes(n);
            
            mesh.nodes(iNode).iBElements = ...
                [mesh.nodes(iNode).iBElements iBElement];
            mesh.nodes(iNode).iBFaces = ...
                [mesh.nodes(iNode).iBFaces iFace];
        end
    end
end
%% END - Capture boundary elements(faces)


%% BEG - Set some boundary node properties. 
%
for iBoundary=1:mesh.numberOfBoundaries
    
    boundary = mesh.boundaries(iBoundary);
    iFaces = boundary.startFace:boundary.endFace;
    
    mesh.boundaries(iBoundary).owners = ...
        [mesh.faces(iFaces).iOwner];
    
    mesh.boundaries(iBoundary).iNodes = ...
        unique([ mesh.faces(iFaces).iNodes ]);
    
    mesh.boundaries(iBoundary).iFaces = iFaces;
    
    mesh.boundaries(iBoundary).pointFaces = ...
        calculatePointFaces(mesh, iBoundary);
    
    mesh.boundaries(iBoundary).pointNormals = ...
        calculatePointNormals(mesh, iBoundary);
end
%% END - Set some boundary node properties.

end
