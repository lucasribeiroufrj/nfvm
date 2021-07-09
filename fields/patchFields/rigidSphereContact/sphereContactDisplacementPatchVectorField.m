classdef sphereContactDisplacementPatchVectorField < patchVectorField
    %sphereContactDisplacementPatchVectorField Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Nodes (vetices) of the mesh;
        nodes_;
        
        % Used for interpolation (from faces' centroid to node)
        weights_;
        
        % The colliding sphere.
        rigidSphereContact_;
        
        % timeIndex
        timeIndex_;
    end
    
    methods
        function obj = ...
                sphereContactDisplacementPatchVectorField...
                (...
                    mesh,...
                    boundary...
                )
            %sphereContactDisplacementPatchVectorField Construct an instance of this class
            %   Detailed explanation goes here
            
            args{1} = mesh;
            args{2} = boundary.index;
            
            obj@patchVectorField(args{:});
            
            obj.rigidSphereContact_ = ...
                rigidSphereContact...
                (...
                    boundary.radii,...
                    boundary.initialPosition,...
                    boundary.initialVelocity,...
                    obj.mesh().timeMesh().endTime()...
                );
            
            obj = obj.initialize();
        end
        
        function obj = initialize(obj)
            
            interval = obj.startFace():obj.endFace();
            obj.nodes_ = unique([obj.mesh().faces(interval).iNodes]);
            nNodes = length(obj.nodes_);
            weights(nNodes) = struct('total',0,'factors',[],'bFaces',[]);
            mesh = obj.mesh().spaceMesh();
            
            for i = 1:length(obj.nodes_)
                
                node = mesh.nodes(obj.nodes_(i));
                nodeCentroid = node.centroid;
                
                % Get the faces that belong to the boundary patch
                faces = ...
                    node.iFaces...
                    (...
                          node.iFaces >= obj.startFace() ...
                        & node.iFaces <= obj.endFace()...
                    );
                
                weights(i).total = 0;
                
                for idx = 1:length(faces)

                    iFace = faces(idx);
                    face = mesh.faces(iFace);
                    c = face.centroid;
                    d_vc = nodeCentroid - c;
                    invDist = 1 /norm(d_vc);
                    weights(i).factors = [weights(i).factors invDist];
                    weights(i).total = weights(i).total + invDist;
                    bFace = iFace - mesh.numberOfInteriorFaces;
                    weights(i).bFaces = [weights(i).bFaces  bFace];
                end
            end
            
            obj.weights_ = weights;
            
            obj.timeIndex_ = -1;
        end
        
        function obj = updateDivImpKgradUfluxes(obj, solidModel)
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            obj = obj.correct();
            obj = obj.updateDivImpKgradUfluxes@patchVectorField(solidModel);
        end
        
        function obj = correct(obj, volU, volGradU, surfaceGradU) %#ok<INUSD>
            %correct Computes the dirichlet (displacement) boundary 
            % condition, i.e:
            %
            %   For each vertex v:
            %       Check if there is a collision (vertex of the mesh is 
            %       inside the sphere).
            %       Find the penetration vector deltaU.
            %       Get the faces connected to vertex.
            %       for each connected face f:
            %           Add this deltaU to the faces's centroid
            %           displacement field (using interpolation).
            %
            %   U_b stands for displacement at the boundary.
            %
            
            currentTimeIndex = obj.mesh().timeMesh.timeIndex();
            if obj.timeIndex_ ~= currentTimeIndex
                
                obj.timeIndex_ = currentTimeIndex;
            else
                return;
            end
            
            mesh = obj.mesh().spaceMesh();
            startFace = obj.startFace();
            endFace = obj.endFace();
            currentTime = obj.mesh().timeMesh().value();
            Y = obj.rigidSphereContact_.position(currentTime);
            sphereCentroid = [ 0; Y; 0 ];
            radii = obj.rigidSphereContact_.radii();
            
            U_b = vectorField(obj.numberOfFaces());
            factor = scalarField(obj.numberOfFaces);
            
            for i = 1:length(obj.nodes_)
                
                node = mesh.nodes(obj.nodes_(i));
                
                % Interpolation for U (from faces' centroid to node). This
                % calculates the current position of node.
%                 dU = 0;
%                 weight = obj.weights_(i);
%                 for j = 1:length(weight.factors)
%                     
%                     U_f = obj.field(:,weight.bFaces(j));
%                     dU = dU + weight.factors(j)*U_f;
%                 end
%                 dU = dU/weight.total;
                
                %nodeCentroid = node.centroid + dU;
                nodeCentroid = node.centroid;
                diffVector = nodeCentroid - sphereCentroid;
                
                % We need only to check norm squared
                dist = norm( diffVector, 2);
                
                if dist < radii
                    % Collided, so compute the penetration
                    
                    len = norm(diffVector);
                    depth = radii - len;
                    direction = diffVector/len;
                    deltaU = depth*direction;
                    
                    %% Add this deltaU contribution to the face's centroid 
                    % displacement.
                    
                    % Get the faces that belong to the boundary patch
                    faces = ...
                        node.iFaces...
                        (...
                              node.iFaces >= startFace ...
                            & node.iFaces <= endFace...
                        );
                    
                    for idx = 1:length(faces)
                        
                        iFace = faces(idx);
                        face = mesh.faces(iFace);
                        c = face.centroid;
                        d_vc = nodeCentroid - c;
                        invDist = 1 /norm(d_vc);
                        l = iFace - startFace + 1;
                        U_b.data(:,l) = U_b.range(l) + invDist*deltaU;
                        factor.data(l) = factor.range(l) + invDist;
                    end
                    %%
                end
            end
            
            % Avoids zero-division.
            for i = 1:factor.end(1,1)
                if factor.range(i) == 0
                    factor.data(i) = 1;
                end
            end
            
            obj = obj.setVectorField(U_b / factor);
        end
        
        function obj = update(obj)
            
            obj = obj.correct();
        end
    end
end

