classdef CantileverBeam2d < Problem & handle
    %CantileverBeam Summary of this class goes here
    %   Detailed explanation goes here
    %
    % Case described in [4.2; Demirdzic - 1994].
    %
    % Ref.:
    %   [Demirdzic - 1994] - Finite volume method for stress analysis in 
    % complex domains
    
    properties
        
        % Gravity vector (body force)
        gravity_ = [ 0; -9.8; 0 ];
        
        % To draw to mesh
        meshDrawer_;
        
        % Boundary conditions for displacement U 
        boundaryConditions_;
    end
    
    methods
       
        function obj = CantileverBeam2d(casePath, params)
            %CantileverBeam2d Construct the case described in 
            %   The beam is fixed on the left.
            
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                    'nTimeSteps', 0 ...
                );
            
            obj@Problem(setup, params);
            
            obj.setBoundaryConditions();
            obj.setMeshDrawer();
            obj.setSolidModel(params);
        end
        
        function material = createMaterial...
                (...
                    obj,...
                    materialEnum,...
                    varargin...
                )
            
            material = ...
                materials.createMaterial...
                (...
                    materialEnum,...
                    planeStress(false),...
                    varargin{:}...
                );
        end
        
        function obj = setBoundaryConditions(obj)
           
            for iBoundary=1:obj.mesh().numberOfBoundaries
    
                boundary = obj.mesh().boundaries(iBoundary);

                bName = boundary.userName;
                boundCondition = [];
                
                if strcmp(bName,'top') || strcmp(bName,'bottom')

                    % Traction free boundary
                    boundCondition.kind = 'solidTraction';
                    
                elseif strcmp(bName,'right')
                    
                    boundCondition.kind = 'solidTraction';
                    
                    nFaces = boundary.numberOfBFaces();
                    data = zeros(3, nFaces);
                    
                    for i=1:nFaces
                        
                        data(:,i) = [ 0; -1e+6; 0; ];
                    end
                    
                    boundCondition.traction = vectorField(data);
                
                elseif strcmp(bName,'left')

                    % The base is fixed on the ground.
                    boundCondition.kind = 'fixedDisplacement';
                    
                elseif strcmp(bName,'front') || strcmp(bName,'back')
                    
                    % Set as 2d case.
                    boundCondition.kind = 'empty';
                else
                    ME = MException('CantileverBeam2d:boundaryCondition',...
                        'Invalid boundary name: %s', bName);
                    throw(ME)
                end
                
                boundCondition.index = iBoundary;
                obj.boundaryConditions_{iBoundary} = boundCondition;
            end
        end
        
        function obj = setMeshDrawer(obj)
            
            obj.meshDrawer_ = BoundaryMeshDrawer(obj.mesh());
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
        
        function meshDrawer = meshDrawer(obj)
            
            meshDrawer = obj.meshDrawer_;
        end
        
        function setCamera(obj)
            
            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],...
                'XLim',[-0.01 2.01],...
                'YLim',[-0.3 0.11],...
                'ZLim',[-0.01 1.01]...
             )
            xlabel('x') 
            ylabel('y')
            zlabel('z')
        end
        
        function bool = hasBodyForce(obj)
           
            %bool = true;
            bool = false;
        end
        
        function RHS = bodyForce(obj)
            
            mesh = obj.mesh();
            N = mesh.numberOfElements();
            
            density = obj.materialLaw().density();
            Density_g = vectorField(N);
            density_g = density * obj.gravity_;
            Density_g(1:end) = density_g;
            
            force = mesh.volume() * Density_g;
            
            RHS = -force.data.';
        end
        
        function maxDisplacementError...
                (...
                    obj,...
                    numericalVolU,...
                    numericalVolGradU...
                )
            
            bIndex = 3; % Right
            solution = -0.01456;
            
            if isa(obj.solidModel(), 'Segregated')
            
                mesh = obj.mesh();

                CN = mesh.CN();
                boundary = CN.boundaryField(bIndex);
                iFace = boundary.startFace();
                iOwner = mesh.faces(iFace).iOwner;
                CN = boundary.field(1);
                U = numericalVolU.internalField().field(iOwner);
                gradU = numericalVolGradU.internalField().field(iOwner);

                U = U + gradU*[CN(1); 0; 0];
                nHeight = U(2);
                
            elseif isa(obj.solidModel(), 'BlockCoupled') ...
                || isa(obj.solidModel(), 'NonLinearBlockCoupled') ...
                || isa(obj.solidModel(), 'BlockCoupledWithNewBoundary')
                
                U = numericalVolU.boundaryField(bIndex).field(1);
                nHeight = U(2);
            else
                error('No valid solidModel.')
            end
            
            diff = abs(solution - nHeight);
            fprintf('U = %d\n', nHeight);
            sError = sprintf('%0.2g%%', diff/abs(solution)*100);
            fprintf('Beam deflection error: %s: \n', sError);
        end
        
        function obj = statistics(obj)
            
            volU = obj.solidModel.volU();
            volGradU = obj.solidModel.volGradU();
            obj.maxDisplacementError(volU, volGradU);
        end
        
        function solidModelHasEvolved(obj)
            %solidModelHasEvolved should be called after the solidModel has
            % evolved. Override this method as you wish.
            
            obj.statistics();
            
            obj.solidModelHasEvolved@Problem();
        end
        
        function bool = includeAcceleration(obj)
            % includeAcceleration Override this function if you need to
            %   remove the temporal term acceleration. Change the sentence
            %   below to "bool = false".
            
            bool = false;
        end
    end
end

