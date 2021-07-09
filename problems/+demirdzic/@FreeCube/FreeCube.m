classdef FreeCube < Problem & handle
    %FreeCube Summary of this class goes here
    %   Detailed explanation goes here
    %
    % Case described in [4.1.a; Demirdzic - 1994].
    %
    % Ref.:
    %   [Demirdzic - 1994] - Finite volume method for stress analysis in 
    % complex domains
    
    properties
        
        % To draw to mesh.
        meshDrawer_;
        
        % Boundary conditions for displacement U .
        boundaryConditions_;
        
        % Variables of the problem.
        alpha_ = 2;     % Random value;
        deltaTemp_ = 3;   % Random value;
    end
    
    methods
        
        function obj = FreeCube(casePath, params)
            %FreeCube Construct the case described in 
            %   [4.1.a; Demirdzic - 1994].
            
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                    'nTimeSteps', 0 ... % Static case
                );
            
            obj@Problem(setup, params);
            
            numberOfCorrectors = ...
                tryGetOrDefault(params, 'numberOfCorrectors', 50);
            
            obj.setBoundaryConditions();
            
            obj.setMeshDrawer();
            
            relaxation = ...
                tryGetOrDefault(params, 'relaxation', 0.90);
            
            obj.setRelaxation(relaxation);
            obj.setNumberOfCorrectors(numberOfCorrectors);
            obj.setSolutionTolerance(1e-9);

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
                
                if strcmp(bName,'bottom') || strcmp(bName,'left') ...
                    || strcmp(bName,'back')
                    
                    boundCondition.kind = 'solidSymmetry';
                    
                elseif strcmp(bName,'top') || strcmp(bName,'right') ...
                    || strcmp(bName,'front')
            
                    % Traction free boundary
                    boundCondition.kind = 'solidTraction';
                else
                    ME = MException('FreeCube:setBoundaryCondition',...
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
        
        function materialLaw = createMaterialLaw...
                (...
                    obj,...
                    materialLawEnum,...
                    material,...
                    varargin...
                )
            
            materialLaw = ...
                materialLaws.createMaterialLaw...
                (...
                    materialLawEnum,...
                    obj.mesh(),...
                    material,...
                    obj.alpha_,...
                    obj.deltaTemp_...
                );
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
            (   gca,...
                'CameraPosition', [-12.4 -16.3 11.4],...
                'CameraUpVector', [0.18 0.23 0.95],...
                'CameraViewAngle', 20.32,...
                'DataAspectRatio', [1 1 1],...
                'PlotBoxAspectRatio', [1 1 1],...
                'CameraTarget', [3.15 3.36 3.57]...
            );
        end
        
        function maxDisplacementError(obj, numericalVolUField)
            
            mesh = obj.mesh();
            maxU = zeros(1,3);
            maxUAnalytical = zeros(1,3);
            
            U = numericalVolUField.internalField();%(:,iElem);
            
            for iElem = 1:mesh.numberOfElements
                
                element = mesh.elements(iElem);
                numerical = U.range(iElem);
                analytical = obj.alpha_*obj.deltaTemp_*element.centroid;
                diffPhi = abs(numerical - analytical);
                
                vec = [ diffPhi(1,1) diffPhi(1,1) diffPhi(2,1) ];
                sta = [ analytical(1,1) analytical(1,1) analytical(2,1) ];
                
                for i = 1:3
                    
                    if vec(i) > maxU(i)
                        
                        maxU(i) = vec(i);
                        maxUAnalytical(i) = abs(sta(i));
                    end
                end
            end
            
            U = sprintf('U = %0.2e (%0.2g%%)', vec(1), ...
                (maxU(1)/(maxUAnalytical(1) + eps))*100);
            W = sprintf('W = %0.2e (%0.2g%%)', vec(2), ...
                (maxU(2)/(maxUAnalytical(2) + eps))*100);
            Z = sprintf('Z = %0.2e (%0.2g%%)', vec(3), ...
                (maxU(3)/(maxUAnalytical(3) + eps))*100);
            
            fprintf('Max disp. error (m):\n%s\n%s\n%s\n', U, W, Z);
        end
        
        function obj = statistics(obj)
            
            volU = obj.solidModel.volU();
            obj.maxDisplacementError(volU);
        end
        
        function bool = includeAcceleration(obj)
            % includeAcceleration Override this function if you need to
            %   remove the temporal term acceleration. Change the sentence
            %   below to "bool = false".
            
            bool = false;
        end
    end
end

