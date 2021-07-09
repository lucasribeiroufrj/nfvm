function mesh = fvmReadOpenFoamMesh(caseDirectory)

    function instance = allocateBoundaryContainer(numberOfBoundaries)
        
        instance = repmat(struct(...
            'userName', '',...
            'index', inf,...
            'type', '',...
            'startFace', inf,...
            'numberOfBFaces', inf,...
            'endFace', inf...
            ),...
            1,...
            numberOfBoundaries);
        
        range = num2cell(1:numberOfBoundaries); 
        [instance.index] = range{:};
    end

    function instance = allocateElementsContainer(numberOfElements)
        
        instance = repmat(struct(...
            'index', inf,...
            'iNeighbours', [],...
            'iFaces', [],...
            'iNodes', [],...
            'volume', 0,...
            'faceSign', []...
            ),...
            1,...
            numberOfElements);
        
        range = num2cell(1:numberOfElements);
        [instance.idex] = range{:};
    end

if(nargin==0)
    caseDirectory = '~';
    caseDirectory = uigetdir(caseDirectory);
end

pointsFile = [caseDirectory,'/constant/polyMesh/points'];
facesFile = [caseDirectory,'/constant/polyMesh/faces'];
ownerFile = [caseDirectory,'/constant/polyMesh/owner'];
neighbourFile = [caseDirectory,'/constant/polyMesh/neighbour'];
boundaryFile = [caseDirectory,'/constant/polyMesh/boundary'];

%
%readPoints
%
fpid = fopen(pointsFile, 'r');
for i=1:18
    fgetl(fpid);
end
numberOfPoints = fscanf(fpid, '%d',1);
fgetl(fpid);
fgetl(fpid);
for n = 1:numberOfPoints
    fscanf(fpid, '%c',1);
    x=fscanf(fpid, '%f',1);
    y=fscanf(fpid, '%f',1);
    z=fscanf(fpid, '%f',1);
    fscanf(fpid, '%c',1);
    fscanf(fpid, '%c',1);
    nodes(n).centroid = [x y z]';
    nodes(n).index = n;
    nodes(n).iFaces=[];
    nodes(n).iElements=[];
end
mesh.nodes = nodes;
mesh.numberOfNodes = numberOfPoints;

mesh.caseDirectory = caseDirectory;
fclose(fpid);

%
%read Faces
%
ffid = fopen(facesFile, 'r');
%header = textscan(fpid,'%s',18,'EndOfLine','\n') 
numberOfFaces = textscan(ffid,'%d',1,'HeaderLines',18);
numberOfFaces = numberOfFaces{1};
header = fscanf(ffid, '%s', 1);

% read each set of measurements
for n = 1:numberOfFaces
    theNumberOfPoints = fscanf(ffid, '%d', 1);
    header = fscanf(ffid, '%c', 1);
    faces(n).iNodes=fscanf(ffid, '%d', theNumberOfPoints)'+1;
    header = fscanf(ffid, '%c', 1);
    faces(n).index=n;
    faces(n).iOwner = -1;
    faces(n).iNeighbour = -1;
end
mesh.numberOfFaces = numberOfFaces;
fclose(ffid);

%
%read owners
%
 foid = fopen(ownerFile, 'r');
%header = textscan(fpid,'%s',18,'EndOfLine','\n') 
numberOfOwners = textscan(foid,'%d',1,'HeaderLines',18);
numberOfOwners = numberOfOwners{1};
header = fscanf(foid, '%s', 1);
for n=1:numberOfOwners
    faces(n).iOwner = fscanf(foid,'%d',1)+1;
    faces(n).iNeighbour=-1;
end
fclose(foid);
%
%read neighbours
%
fnid = fopen(neighbourFile, 'r');
%header = textscan(fpid,'%s',18,'EndOfLine','\n') 
numberOfInteriorFaces = textscan(fnid,'%d',1,'HeaderLines',18);
numberOfInteriorFaces = numberOfInteriorFaces{1};
header = fscanf(fnid, '%s', 1);
for n=1:numberOfInteriorFaces
    faces(n).iNeighbour = fscanf(fnid,'%d',1)+1;
end

mesh.faces = faces;
mesh.numberOfInteriorFaces = numberOfInteriorFaces;
mesh.numberOfElements = max([faces.iNeighbour]);
fclose(fnid);

%
%read boundaries
%
fbid = fopen(boundaryFile, 'r');
numberOfBoundaries = textscan(fbid,'%d',1,'HeaderLines',17);
numberOfBoundaries = numberOfBoundaries{1};
fscanf(fbid, '%c', 2); % Ignore the header

boundaries = allocateBoundaryContainer(numberOfBoundaries);

for n=1:numberOfBoundaries
    
    boundaries(n).userName = fscanf(fbid, '%s', 1);
    token='{';
    
    while(strcmp(token,'}')==false)
        
        token = fscanf(fbid, '%s', 1);
        
        if(strcmp(token,'type'))
            
            theType = fscanf(fbid, '%s', 1);
            boundaries(n).type = theType(1:length(theType)-1);   
        end
        
        if(strcmp(token,'startFace'))
            
            boundaries(n).startFace = fscanf(fbid, '%d', 1)+1;   
        end
        
        if(strcmp(token,'nFaces'))
            
            boundaries(n).numberOfBFaces = fscanf(fbid, '%d', 1);   
        end
        
        fgetl(fbid);
    end
    
    boundaries(n).endFace = ...
        boundaries(n).startFace + boundaries(n).numberOfBFaces - 1;
end
mesh.boundaries = boundaries;
mesh.numberOfBoundaries = numberOfBoundaries;
mesh.numberOfPatches = numberOfBoundaries;
fclose(fbid);

elements = allocateElementsContainer(mesh.numberOfElements);

%% BEG - Initialize Finite Volume Elements 
%
for iFace=1:mesh.numberOfInteriorFaces
    
   iOwner = faces(iFace).iOwner;
   iNeighbour = faces(iFace).iNeighbour;

   elements(iOwner).iFaces = [elements(iOwner).iFaces iFace];
   elements(iOwner).faceSign = [elements(iOwner).faceSign 1];
   elements(iOwner).iNeighbours = [elements(iOwner).iNeighbours iNeighbour];

   elements(iNeighbour).iFaces = [elements(iNeighbour).iFaces iFace];
   elements(iNeighbour).faceSign = [elements(iNeighbour).faceSign -1];
   elements(iNeighbour).iNeighbours = [elements(iNeighbour).iNeighbours iOwner];
end

for iBFace=mesh.numberOfInteriorFaces+1:mesh.numberOfFaces
    
   iOwner = mesh.faces(iBFace).iOwner;

   elements(iOwner).iFaces= [elements(iOwner).iFaces iBFace];
   elements(iOwner).faceSign =[elements(iOwner).faceSign 1];
end

for iElement=1:mesh.numberOfElements
    
   elements(iElement).numberOfNeighbours = length(elements(iElement).iNeighbours);
end
%% END - Initialize elements

numberOfFaces = mesh.numberOfFaces;
mesh.elements = elements;
mesh.numberOfBElements = numberOfFaces - numberOfInteriorFaces;
mesh.numberOfBFaces = numberOfFaces - numberOfInteriorFaces;
mesh.numberOfTotalElements = mesh.numberOfElements + mesh.numberOfBElements;

mesh=fvmProcessOpenFoamMesh(mesh);

end
