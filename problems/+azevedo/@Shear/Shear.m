classdef Shear < Problem & handle
    %Shear Tangential displacement is applied in on of the boundaries.
    %
    % TODO:
    %   -/-
    %
    % Nomenclature.: 
    % - The letters U stand for displacement vector
    % respectively.
    % - The U_o, U_oo stands for U at old, and old old time respectively.

    
    properties
        
        % The amount of shear.
        deltaShear_
        
        % To draw to mesh
        meshDrawer_;
        
        % Boundary conditions for displacement U 
        boundaryConditions_;
    end
    
    methods
        
        function obj = Shear(casePath, params)
            % Shear Creates a shear unit cube case.
            %   Shear(casePath) Creates a problem which uses the mesh
            %   stored at casePath, a default mechanical law (BlatzKo) and
            %   also a default amount of shearing given by deltaLength 
            %   (0.45).
            
            nTimeSteps = ...
                tryGetOrDefault(params, 'numberOfTimeSteps', 1);
            
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                    'nTimeSteps', nTimeSteps,...
                    'startTime', 0.0,...
                    'endTime', 1.0...
                );
            
            obj@Problem(setup, params);
            
            obj.deltaShear_ = tryGetOrDefault(params, 'deltaLength', 0.45);
            
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
                
                if strcmp(bName,'right') || strcmp(bName,'left') ...
                        || strcmp(bName,'front') || strcmp(bName,'back')

                    % Traction free boundary
                    boundCondition.kind = 'solidTraction';

                elseif strcmp(bName,'top')

                    boundCondition.kind = 'uniformVariableDisplacement';
                    
                    boundCondition.displacementFunction = ...
                        @(time) obj.deltaLength(time);
                
                elseif strcmp(bName,'bottom')

                    % The base is fixed on the ground.
                    boundCondition.kind = 'fixedDisplacement';
                else
                    ME = MException('Shear:boundaryCondition',...
                        'Invalid boundary name: %s', bName);
                    throw(ME)
                end
                
                boundCondition.index = iBoundary;
                obj.boundaryConditions_{iBoundary} = boundCondition;
            end
        end
        
        function value = deltaLength(obj, time)
            
            nf = obj.mesh().nf().boundaryField(3);
            value = (time*obj.deltaShear_)*nf;
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
            
            xLim = 1.1;
            
            if obj.deltaShear_ > 0
                xLim = 1 + obj.deltaShear_ + 0.1;
            end
            
            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],...
                'XLim',[-0.1 xLim],...
                'YLim',[-0.1 1.1],...
                'ZLim',[-0.1 1.1]...
             )
            xlabel('x') 
            ylabel('y')
            zlabel('z')
            
            camorbit(20, 20,'data',[0 1 0]);
        end
    end
end

