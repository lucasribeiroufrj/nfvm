classdef fvMesh < handle
    %fvMesh FVM Space-time mesh handle 
    %   Detailed explanation goes here
    %
    % TODO:
    % 
    
    properties
        
        % Time handle.
        timeMesh_;
        
        % Handle to the geometric mesh.
        spaceMesh_;
        
        % Volume identity matrix field.
        volI_;
        
        % Surface identity matrix field.
        surfaceI_;
        
        % Three unit base vectors (e1,e2,e3) as volume field.
        volUnitBaseVectors_;
        
        % Three unit base vectors (e1,e2,e3) as surface field.
        surfaceUnitBaseVectors_;
        
        % Surface vector field.
        Sf_;
        
        % Surface area = |Sf_| field.
        SfNorm_;
        
        % The a surface field Sf_ / SfNorm_.
        nf_;
        
        % It stores each distance vector from the owner element centroid to 
        % the surface centroid. It is a surface field.
        CN_;
        
        % The surface field |CN_|.
        CNnorm_;
        
        % The surface field CN_ / |CN_|.
        eCN_;
        
        % Volume centroid field (for internal field plus face centroid 
        % boundary fields)
        r_C_;
        
        % A surface field storing the face centroid.
        r_Cf_;
        
        % Stores the volume of each finite volume (as a scalarField)
        volume_;
        
        % Face Weighting Factors (or Geometric factors).
        % See [pg. 159, Moukalled - 2015].
        wf_;
        wc_;
        
        % Vector joining the centroids of the elements stradding the
        % interface. This is surfaceVectorField.
        dCF_ ;
        
        % |dCF_|
        dCFnorm_;
        
        % dCF_ / dCFnorm_
        eCF_;
        
        % Boudnary information.
        boundaries_;
        
        % One-to-one face-owner and face-neighbour lists (internal faces 
        % only). This allows parallel assignment.
        faceOwnerLists_;
        faceNeighLists_;
        
        % The iOwners_(i) == cell that owns face i.
        % The iNeighs_(i) == cell that neighbours face i.
        % This is only for interior faces.
        iOwners_;
        iNeighs_;
    end
    
    methods
        function obj = fvMesh(timeMesh, spaceMesh)
            %Mesh Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.timeMesh_ = timeMesh;
            obj.spaceMesh_ = spaceMesh;
            
            obj.processSpaceMesh();
        end
        
        %% BEG - Data member access
        function timeMesh = timeMesh(obj)
            
            timeMesh = obj.timeMesh_;
        end
        
        function spaceMesh = spaceMesh(obj)
            
            spaceMesh = obj.spaceMesh_;
        end
        
        function obj = setSpaceMesh(obj, spaceMesh)
            
            obj.spaceMesh_ = spaceMesh;
        end
        
        function volI = volI(obj)
            
            volI = obj.volI_;
        end
        
        function surfaceI = surfaceI(obj)
            
            surfaceI = obj.surfaceI_;
        end
        
        function value = volUnitBaseVectors(obj)
            
            value = obj.volUnitBaseVectors_;
        end
        
        function value = surfaceUnitBaseVectors(obj)
            
            value = obj.surfaceUnitBaseVectors_;
        end
        
        function Sf = Sf(obj)
            
            Sf = obj.Sf_;
        end
        
        function SfNorm = SfNorm(obj)
            
            SfNorm = obj.SfNorm_;
        end
        
        function nf = nf(obj)
            
            nf = obj.nf_;
        end
        
        function CN = CN(obj)
            
            CN = obj.CN_;
        end
        
        function eCN = eCN(obj)
            
            eCN = obj.eCN_;
        end
        
        function CNnorm = CNnorm(obj)
            
            CNnorm = obj.CNnorm_;
        end
        
        function r_C = r_C(obj)
            r_C = obj.r_C_;
        end
        
        function r_Cf = r_Cf(obj)
            r_Cf = obj.r_Cf_;
        end
        
        function volume = volume(obj)
            
            volume = obj.volume_;
        end
        
        function wf = wf(obj)
            
            wf = obj.wf_;
        end
        
        function wc = wc(obj)
            
            wc = obj.wc_;
        end
        
        function dCFnorm = dCFnorm(obj)
           
            dCFnorm = obj.dCFnorm_;
        end
        
        function  eCF = eCF(obj)
            
            eCF = obj.eCF_;
        end
        
        function boundary = boundaries(obj, index)
           
            boundary = obj.boundaries_{index};
        end
        
        function faceOwnerLists = faceOwnerLists(obj)
           
            faceOwnerLists = obj.faceOwnerLists_;
        end
        
        function faceNeighLists = faceNeighLists(obj)
           
            faceNeighLists = obj.faceNeighLists_;
        end
        
        function iOwners = iOwners(obj)
           
            iOwners = obj.iOwners_;
        end
        
        function iNeighs = iNeighs(obj)
           
            iNeighs = obj.iNeighs_;
        end
        %% END 
        
        function boundary = setBoundary(obj, index, boundary)
            
            obj.boundaries_{index} = boundary;
        end
        
        function nodes = nodes(obj, index)
           
            nodes = obj.spaceMesh_.nodes(index);
        end
        
        function numberOfNodes = numberOfNodes(obj)
           
            numberOfNodes = obj.spaceMesh_.numberOfNodes;
        end
        
        function caseDirectory = caseDirectory(obj)
           
            caseDirectory = obj.spaceMesh_.caseDirectory;
        end
        
        function numberOfFaces = numberOfFaces(obj)
           
            numberOfFaces = obj.spaceMesh_.numberOfFaces;
        end
        
        function numberOfElements = numberOfElements(obj)
           
            numberOfElements = obj.spaceMesh_.numberOfElements;
        end
        
        function numberOfTotalElements = numberOfTotalElements(obj)
            
            numberOfTotalElements = obj.spaceMesh_.numberOfTotalElements;
        end
        
        function faces = faces(obj, index)
           
            faces = obj.spaceMesh_.faces(index);
        end
        
        function numberOfInteriorFaces = numberOfInteriorFaces(obj)
           
            numberOfInteriorFaces = obj.spaceMesh_.numberOfInteriorFaces;
        end
        
        function numberOfBoundaries = numberOfBoundaries(obj)
           
            numberOfBoundaries = obj.spaceMesh_.numberOfBoundaries;
        end
        
        function elements = elements(obj, index)
           
            elements = obj.spaceMesh_.elements(index);
        end
        
        function numberOfBElements = numberOfBElements(obj)
           
            numberOfBElements = obj.spaceMesh_.numberOfBElements;
        end
        
        function numberOfBFaces = numberOfBFaces(obj)
           
            numberOfBFaces = obj.spaceMesh_.numberOfBFaces;
        end
        
        function obj = processSpaceMesh(obj)
            
            mesh = obj.spaceMesh_;
            
            %% Adding boundary data
            for b=1:mesh.numberOfBoundaries

                boundary = mesh.boundaries(b);
                startFace = boundary.startFace;
                nBoundaryFaces = boundary.numberOfBFaces;
                endFace = boundary.endFace;
                indices = startFace:endFace;
                                
                %% Adding primitive data
                obj.boundaries_{b}.userName = boundary.userName;
                obj.boundaries_{b}.index = boundary.index;
                obj.boundaries_{b}.type = boundary.type;
                obj.boundaries_{b}.numberOfBFaces = nBoundaryFaces;
                obj.boundaries_{b}.startFace = startFace;
                
                %% Adding extended data
                obj.boundaries_{b}.endFace = endFace;
                obj.boundaries_{b}.owners = [mesh.faces(indices).iOwner];
            end
            %%
            
            %% BEG - Adding extended data
            
            % beg - Create identity matrices
            nTotalElements = mesh.numberOfTotalElements;
            I = eye(3,3);
            rawTensors = zeros(3,3,nTotalElements);
            for i = 1:nTotalElements
               rawTensors(:,:,i) = I; 
            end
            rawTensorField = tensorField(rawTensors);
            obj.volI_ = volTensorField(obj, rawTensorField);
            
            nFaces = mesh.numberOfFaces;
            rawTensors = zeros(3,3,nFaces);
            for i = 1:nFaces
               rawTensors(:,:,i) = I; 
            end

            rawTensorField = tensorField(rawTensors);
            obj.surfaceI_ = surfaceTensorField(obj, rawTensorField);
            % end
            
            
            % Beg - Create volUnitBaseVectors.
            nTotalElements = mesh.numberOfTotalElements;
            nTotalFaces = mesh.numberOfFaces;
            for i = 1:3
                
                rawVectors = zeros([3,nTotalElements]);
                rawVectors(i,:) = 1;
                obj.volUnitBaseVectors_{i} = ...
                    volVectorField(obj, vectorField(rawVectors));
                
                rawVectors = zeros([3,nTotalFaces]);
                rawVectors(i,:) = 1;
                obj.surfaceUnitBaseVectors_{i} = ...
                    surfaceVectorField(obj, vectorField(rawVectors));
            end
            % End - Create volUnitBaseVectors.
            
            
            obj.Sf_ = ...
                surfaceVectorField...
                (...
                    obj,...
                    vectorField([mesh.faces(:).Sf])...
                );

            obj.SfNorm_ = obj.Sf_.norm();

            obj.nf_ = surfaceVectorField...
                (...
                    obj,...
                    vectorField([mesh.faces(:).nf])...
                );
            
            CN = vectorField([mesh.faces(:).CN]);
            obj.CN_ = surfaceVectorField(obj, CN);
            obj.CNnorm_ = obj.CN_.norm();
            
            eCN = vectorField([mesh.faces(:).eCN]);
            obj.eCN_ = surfaceVectorField(obj, eCN);

            centroids = [mesh.elements(:).centroid];
            for b=1:mesh.numberOfBoundaries
               
                boundary = obj.boundaries(b);
                indices = boundary.startFace:boundary.endFace;
                r_Cf = [mesh.faces(indices).centroid ];
                centroids = [ centroids r_Cf ];
            end
            r_C = vectorField(centroids);
            obj.r_C_ = volVectorField(obj, r_C);
            
            obj.volume_ = scalarField([mesh.elements(:).volume]);
            
            r_Cf = vectorField([mesh.faces(:).centroid]);
            obj.r_Cf_ = surfaceVectorField(obj, r_Cf);
            
            %% beg 
            nFaces = mesh.numberOfFaces;
            dCFnorm = zeros(1, nFaces);
            dCF = zeros(3, nFaces);
            eCF = zeros(3, nFaces);
            
            for iFace=1:mesh.numberOfInteriorFaces
                
                face = mesh.faces(iFace);
                iOwner = face.iOwner;
                iNeigh = face.iNeighbour;
                rC = mesh.elements(iOwner).centroid;
                rF = mesh.elements(iNeigh).centroid;
                dcf = rF - rC;
                dcfNorm = norm(dcf);
                dCFnorm(iFace) = dcfNorm;
                dCF(:,iFace) = dcf;
                eCF(:,iFace) = dcf/dcfNorm;
            end
            
            startFace = mesh.numberOfInteriorFaces + 1;
            endFace = startFace + mesh.numberOfBFaces - 1; %
            for iFace = startFace:endFace

                face = mesh.faces(iFace);
                iOwner = face.iOwner;
                %iNeigh = face.iNeighbour; %
                rC = mesh.elements(iOwner).centroid; %
                rF = face.centroid; %
                dcf = rF - rC;
                dcfNorm = norm(dcf);
                dCFnorm(iFace) = dcfNorm;
                dCF(:,iFace) = dcf;
                eCF(:,iFace) = dcf/dcfNorm;
            end
            
            obj.dCFnorm_ = surfaceScalarField(obj, scalarField(dCFnorm));
            obj.dCF_ = surfaceVectorField(obj, vectorField(dCF));
            obj.eCF_ = surfaceVectorField(obj, vectorField(eCF));
            %% end
            
            %% beg 
            WF = [mesh.faces(1:mesh.numberOfFaces).gf];
            wf = scalarField(WF);
            wc = scalarField(1 - WF);
            obj.wf_ = surfaceScalarField(obj, wf);
            obj.wc_ = surfaceScalarField(obj, wc);
            %% end
            
            obj.faceOwnerLists_ = ...
                fvMesh.createOneToOneMapping(mesh, 'iOwner');
            obj.faceNeighLists_ = ...
                fvMesh.createOneToOneMapping(mesh, 'iNeighbour');
            
            iFaces = 1:mesh.numberOfInteriorFaces;
            obj.iOwners_ = [mesh.faces(iFaces).iOwner];
            obj.iNeighs_ = [mesh.faces(iFaces).iNeighbour];
            %% END
        end
        
        function e1 = volE1(obj)
            %e1 Unit base vectors [1 0 0] as volume field.
            
            e1 = obj.volUnitBaseVectors_{1}; 
        end
        
        function e2 = volE2(obj)
            %e2 Unit base vectors [0 1 0] as volume field.
            
            e2 = obj.volUnitBaseVectors_{2}; 
        end
        
        function e3= volE3(obj)
            %e3 Unit base vectors [0 0 1] as volume field.
            
            e3 = obj.volUnitBaseVectors_{3}; 
        end
        
        function e1 = surfaceE1(obj)
            %e1 Unit base vectors [1 0 0] as surface field.
            
            e1 = obj.surfaceUnitBaseVectors_{1}; 
        end
        
        function e2 = surfaceE2(obj)
            %e2 Unit base vectors [0 1 0] as surface field.
            
            e2 = obj.surfaceUnitBaseVectors_{2}; 
        end
        
        function e3= surfaceE3(obj)
            %e3 Unit base vectors [0 0 1] as surface field.
            
            e3 = obj.surfaceUnitBaseVectors_{3}; 
        end
    end
    
    methods(Static)
        
        function lists = createOneToOneMapping(mesh, other)
            % createOneToOneMapping - Creates one-to-one 
            %   face-owner/neighbour lists.
            %   createOneToOneMapping(mesh, 'iOwner'/'iNeigh') Returns a
            %   cell list, say L. Each element has two arrays members,
            %   faces and iOwner, both of the same size. Then the element
            %   ofindex L.iOwner(i) is the owner of the face L.faces(i).
            
            %% BEG - Creates one-to-one face-owner/neighbour lists.
            nFaces = mesh.numberOfInteriorFaces();
            nElements = mesh.numberOfElements();
            tFaces = nFaces;
            faces = 1:nFaces;
            nLists = 1;

            while true

                F = zeros(1, tFaces);
                O = zeros(1, tFaces);
                elements = zeros(1, nElements);

                nCapturedFaces = 0;
                nMissed = 0;

                for idx = 1:tFaces

                    iFace = faces(idx);
                    iElem = mesh.faces(iFace).(other);

                    if elements(iElem) == 0
                        nCapturedFaces = nCapturedFaces + 1;
                        elements(iElem) = -1;
                        F(nCapturedFaces) = iFace;
                        O(nCapturedFaces) = iElem;
                    else
                        nMissed = nMissed + 1;
                        faces(nMissed) = iFace;
                    end

                    if tFaces == idx
                        break
                    end
                end

                lists{nLists}.faces = F(1:nCapturedFaces);
                lists{nLists}.(other) = O(1:nCapturedFaces);

                if nMissed == 0
                    break;
                end

                tFaces = nMissed;
                nLists = nLists + 1;
            end
        end
    end
end

