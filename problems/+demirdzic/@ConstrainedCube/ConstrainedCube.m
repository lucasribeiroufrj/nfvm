classdef ConstrainedCube < Problem & handle
    %ConstrainedCube Summary of this class goes here
    %   Detailed explanation goes here
    %
    % Case described in [4.1.b; Demirdzic - 1994].
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
        
        function obj = ConstrainedCube(casePath, params)
            %ConstrainedCube Construct the case described in 
            %   [4.1.b; Demirdzic - 1994].
            
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                    'nTimeSteps', 0 ... % Static case
                );
            
            obj@Problem(setup, params);
            
            numberOfCorrectors = ...
                tryGetOrDefault(params, 'numberOfCorrectors', 2);
            
            obj.setBoundaryConditions();
            
            obj.setMeshDrawer();
            
            %obj.setRelaxation(0.98);
            obj.setNumberOfCorrectors(numberOfCorrectors);
            obj.setSolutionTolerance(4);
            
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
                    || strcmp(bName,'back') || strcmp(bName,'top') ...
                    || strcmp(bName,'right') || strcmp(bName,'front')
                    
                    % Zero displacement.
                    boundCondition.kind = 'fixedDisplacement';
                else
                    ME = MException('ConstrainedCube:setBoundaryCondition',...
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
                'CameraPosition', [-13.5 -16.5 4.6],...
                'CameraUpVector', [0.09 0.12 0.98],...
                'CameraViewAngle', 6.608,...
                'DataAspectRatio', [1 1 1],...
                'PlotBoxAspectRatio', [1 1 1]...
            );
        end
        
        function maxStressError...
                 (...
                    obj,...
                    numericalVolStress,...
                    analyticalVolStress...
                 )
            
            mesh = obj.mesh();
            
            numericalStress = numericalVolStress.internalField();
            analyticalStress = analyticalVolStress.internalField();
            diff = numericalStress - analyticalStress;
            
            max = 0;
            
            for iElem = 1:mesh.numberOfElements
                
                value = norm(diff.range(iElem), Inf);
                
                if value > max
                    max = value;
                end
            end
            
            sErrorStress = sprintf(['MAX_{for all elements}', ...
                '(|| Analyt.Stress - Numerical.Stress ||_inf) = ', ...
                '%0.2e'], max);
            
            fprintf('Max stress error (Pa): %s\n', sErrorStress);
        end
        
        function obj = statistics(obj)
            
            zeroVolGradU = volTensorField(obj.mesh());
            mModel = obj.solidModel.mechanicalModel();
            analyticalVolStress = mModel.computePiolaField(zeroVolGradU);
            
            volGradU = obj.solidModel.volGradU();
            
            numericalVolStress = mModel.computePiolaField(volGradU);
            
            obj.maxStressError(numericalVolStress, analyticalVolStress);
        end
        
        function bool = includeAcceleration(obj)
            % includeAcceleration Override this function if you need to
            %   remove the temporal term acceleration. Change the sentence
            %   below to "bool = false".
            
            bool = false;
        end
    end
end

