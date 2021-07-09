classdef ImpactOnBrick < Problem & handle
    %ImpactOnBrick A rigid sphere collides with the top boundary.
    %   of brick made of cork [Fernandes - 2014].
    %   The brick has the dimensions: 129.69x122x34.23 (mm). The sphere's 
    % diameter is 94mm.
    %
    % TODO:
    %   -/-
    %
    % Nomenclature.: 
    % - The letters U and V stand for displacement and velocity vectors
    % respectively.
    % - The U_o, U_oo stands for U at old, and old old time respectively.
    % Same applies for V.
    % - P stands for first piola.
    % - S stands for surface area vector.
    % - U_b stands for displacement at the boundary.
    %
    % Ref.:
    % [Fernandes - 2014]  Modelling impact response of agglomerated.
    %
    
    properties
        
        % Impactor
        impactorSphere_;
        
        % To use symmetric broundary contition.
        isSymmetricCase_;
        
        % To draw to mesh
        meshDrawer_;
        
        % Boundary conditions for displacement U 
        boundaryConditions_;
    end
    
    methods
        
        function obj = ImpactOnBrick(casePath, params)
            %ImpactOnBrick Construct rigid-sphere colliding case.
            %   Detailed explanation goes here
            
            t0 = 0.0;
            tf = 9e-3; % 9 ms
            deltaT = 2e-5;
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                    'nTimeSteps', (tf-t0)/deltaT,...
                    'startTime', t0,...
                    'endTime', tf... 
                );
            
            obj@Problem(setup, params);
            
            obj.setIsSymmetricCase();
               
            obj.setBoundaryConditions();
            
            obj.setMeshDrawer();
            
            obj.setAmgRelaxation(0.98);
            
            obj.setNumberOfCorrectors(...
                tryGetOrDefault(params, 'numberOfCorrectors', 150));
            
            obj.setSolidModel(params);
        end
        
        function obj = setBoundaryConditions(obj)
            
            for iBoundary=1:obj.mesh().numberOfBoundaries
    
                boundary = obj.mesh().boundaries(iBoundary);

                if strcmp(boundary.type,'empty')
                    continue
                end

                bName = boundary.userName;
                boundCondition = [];
                
                if (strcmp(bName,'left') || strcmp(bName,'right') ...
                    || strcmp(bName,'front') || strcmp(bName,'back')) ...
                    && ~obj.isSymmetricCase_

                    % Traction free boundary
                    boundCondition.kind = 'solidTraction';

                elseif (strcmp(bName,'right') ...
                        || strcmp(bName,'front')) && obj.isSymmetricCase_

                    % Traction free boundary
                    boundCondition.kind = 'solidTraction';

                elseif (strcmp(bName,'left') || strcmp(bName,'back')) ...
                        && obj.isSymmetricCase_

                    boundCondition.kind = 'solidSymmetry';

                elseif strcmp(bName,'top')

                    boundCondition.kind = 'sphereContactDisplacement';
                    
                    % 96mm diameter converted to meters
                    boundCondition.radii = (94/2) * 1e-3;
                    
                     % mm to meters
                    heightBrick = 34.23 * 1e-3;
                    
                    % Sphere's centroid position.
                    boundCondition.initialPosition = ...
                        boundCondition.radii + heightBrick; 
                    
                    boundCondition.initialVelocity = -2.5;
                    
                elseif strcmp(bName,'bottom')

                    % The base is fixed on the ground.
                    boundCondition.kind = 'fixedDisplacement';
                    
                else
                    ME = MException('nFVM:boundaryCondition',...
                        'Invalid boundary name: %s', bName);
                    throw(ME)
                end
                
                boundCondition.index = iBoundary;
                obj.boundaryConditions_{iBoundary} = boundCondition;
            end
        end
        
        function setIsSymmetricCase(obj)
            
            [~,name] = ...
                fileparts(obj.mesh().spaceMesh().caseDirectory());
            
            if isempty(strfind(name,'Sym'))
                
                obj.isSymmetricCase_ = false;
            else
                obj.isSymmetricCase_ = true;
            end
        end
        
        function obj = setMeshDrawer(obj)
           
            if obj.isSymmetricCase_
                
                obj.meshDrawer_ = BoundaryMeshDrawer(obj.mesh());
            else
                obj.meshDrawer_ = DrawMeshFacesFirstQuadrant(obj.mesh());
            end
        end
        
        function U = U(obj)
            %U Displacement field
            
            U = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function U_o = U_o(obj)
            %Uo Displacement field at time "old"
            
            U_o = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function U_oo = U_oo(obj)
            %Uoo Displacement field at time "old-old"
            
            U_oo = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function U_b = getDisplacementFromCollision(obj, boundary)
            % getDisplacementFromCollision Computes the dirichlet 
            % (displacement) boundary condition, i.e:
            %
            %   For each vertex v:
            %       Check if there is a collision (vertex is inside the 
            %       sphere).
            %       Find the penetration vector deltaU.
            %       Get the faces connected to vertex.
            %       for each connected face f:
            %           Add this deltaU to the faces's centroid
            %           displacement field (using interpolation).
            %
            %   U_b stands for displacement at the voundary.
            
            theMesh = obj.theMesh();
            startFace = boundary.startFace;
            endFace = startFace + boundary.numberOfBFaces - 1;
            
            currentTime = obj.runTime().value();
            Y = obj.impactorSphere_.position(currentTime);
            sphereCentroid = [ 0; Y; 0 ];
            radii = obj.impactorSphere_.radii();
            % We need to do it only ONCE!
            nodes = unique([theMesh.faces(startFace:endFace).iNodes]);
            
            U_b = vectorField(boundary.numberOfBFaces);
            factor = scalarField(boundary.numberOfBFaces);
            
            for i = 1:length(nodes)
                
                node = theMesh.nodes(nodes(i));
                nodeControid = node.centroid;
                diffVector = nodeControid - sphereCentroid;
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
                        face = theMesh.faces(iFace);
                        c = face.centroid;
                        d_vc = nodeControid - c;
                        invDist = 1 /norm(d_vc);
                        l = iFace - startFace + 1;
                        U_b(l) = U_b(l) + invDist*deltaU;
                        factor(l) = factor(l) + invDist;
                    end
                    %%
                end
            end
            
            % Voids zero-division.
            for i = 1:factor.end(1,1)
                if factor(i) == 0
                    factor(i) = 1;
                end
            end
            
            U_b = U_b / factor;
        end
        
        function meshDrawer = meshDrawer(obj)
            
            meshDrawer = obj.meshDrawer_;
        end
        
        function setCamera(obj)
            
            set...
            (   gca,...
                'CameraPosition', [-0.404 0.273 -0.404],...
                'CameraUpVector', [-0.174 0.974 -0.174],...
                'CameraViewAngle', 6.608,...
                'DataAspectRatio', [1 1 1],...
                'PlotBoxAspectRatio', [1 1 1]...
            );
        end
    end
end

