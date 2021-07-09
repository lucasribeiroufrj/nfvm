classdef Uniaxial < Problem & handle
    %Uniaxial compression/elongation
    %
    % TODO:
    %   - Elogation
    %
    % Nomenclature.: 
    % - The letters U stand for displacement vector
    % respectively.
    % - The U_o, U_oo stands for U at old, and old old time respectively.

    
    properties
        
        % The final length of the cube = original length + deltaLength_.
        deltaLength_
        
        % To draw to mesh
        meshDrawer_;
        
        % Boundary conditions for displacement U 
        boundaryConditions_;
    end
    
    methods
        
        function obj = Uniaxial(casePath, params)
            
            deltaLength = tryGetOrDefault(params, 'deltaLength', 0.45);
            
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
            
            obj.setBoundaryConditions();
             
            obj.setMeshDrawer();
            
            obj.deltaLength_ = deltaLength;
            
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
                
                    
                if strcmp(bName,'top') || strcmp(bName,'bottom') ...
                        || strcmp(bName,'front') || strcmp(bName,'back')

                    % Traction free boundary
                    boundCondition.kind = 'solidTraction';

                elseif strcmp(bName,'right')

                    boundCondition.kind = 'uniformVariableDisplacement';
                    
                    boundCondition.displacementFunction = ...
                        @(time) obj.deltaLength(time);
                
                elseif strcmp(bName,'left')

                    % The base is fixed on the ground.
                    %boundCondition.kind = 'fixedDisplacement';
                    
                    boundCondition.kind = 'solidSymmetry';
                else
                    ME = MException('Uniaxial:boundaryCondition',...
                        'Invalid boundary name: %s', bName);
                    throw(ME)
                end
                
                boundCondition.index = iBoundary;
                obj.boundaryConditions_{iBoundary} = boundCondition;
            end
        end
        
        function value = deltaLength(obj, time)
            
%             nf = obj.mesh().nf().boundaryField(3);
% 
%             value = vectorField((time*obj.deltaLength_)*nf);

            mesh = obj.mesh();
            nf = mesh.nf().boundaryField(3);
            solidModel = obj.solidModel();
            
            % Is incremental approach?
            if false && strfind(lower(class(solidModel)), ...
                    lower('BlockCoupledWithNewBoundary'))
                
                if mesh.timeMesh.isFirstTimeStep()
                    t_o = mesh.timeMesh.startTime();
                else
                    t_o = mesh.timeMesh.valueOld();
                end
                t_ = mesh.timeMesh.value();
                
                increment = (t_ - t_o)*obj.deltaLength_;
            else
                increment = time*obj.deltaLength_;
            end
            
            value = increment*nf;
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
            
            xLim = 1.2;
            
            if obj.deltaLength_ > 0
                xLim = 1 + obj.deltaLength_ + 0.1;
            end
            
            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],...
                'XLim',[-0.2 xLim],...
                'YLim',[-0.2 1.2],...
                'ZLim',[-0.2 1.2]...
             )
            xlabel('x') 
            ylabel('y')
            zlabel('z')
            
            camorbit(20, 20,'data',[0 1 0]);
        end
    end
end

