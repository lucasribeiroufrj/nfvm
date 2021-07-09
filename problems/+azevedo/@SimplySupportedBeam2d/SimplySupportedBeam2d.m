classdef SimplySupportedBeam2d < Problem & handle
    %CantileverBeam Summary of this class goes here
    %   Detailed explanation goes here
    %
    % Case described in [4.2; Demirdzic - 1994].
    %
    % Ref.:
    %   [Demirdzic - 1994] - Finite volume method for stress analysis in 
    % complex domains
    
    properties
        
        % Store Hookean material law.
        materialLaw_;
        
        % Gravity vector (body force)
        gravity_ = [ 0; -9.8; 0 ];
        
        % To draw to mesh
        meshDrawer_;
        
        % Boundary conditions for displacement U 
        boundaryConditions_;
    end
    
    methods
       
        function obj = SimplySupportedBeam2d(casePath, params)
            %SimplySupportedBeam2d Construct the case described in 
            %   The beam is fixed on the left.
            
            nTimeSteps = ...
                tryGetOrDefault(params, 'numberOfTimeSteps', 5);
            
            endTime = ...
                tryGetOrDefault(params, 'endTime', 1.0);
            
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                     'nTimeSteps', nTimeSteps,...
                    'startTime', 0.0,...
                    'endTime', endTime...
                );
            
            obj@Problem(setup, params);
            
            obj.setBoundaryConditions();
            
            obj.setMeshDrawer();
            
            obj.setSolidModel(params);
        end
        
        function obj = setBoundaryConditions(obj)
           
            for iBoundary=1:obj.mesh().numberOfBoundaries
    
                boundary = obj.mesh().boundaries(iBoundary);

                bName = boundary.userName;
                boundCondition = [];
                
                if strcmp(bName,'top') || strcmp(bName,'bottom') ...
                        || strcmp(bName,'right')

                    % Traction free boundary
                    boundCondition.kind = 'solidTraction';
                
                elseif strcmp(bName,'left')

                    % The base is fixed on the ground.
                    boundCondition.kind = 'fixedDisplacement';
                    
                elseif strcmp(bName,'front') || strcmp(bName,'back')
                    
                    % Set as 2d case.
                    boundCondition.kind = 'empty';
                else
                    ME = MException('SimplySupportedBeam2d:boundaryCondition',...
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
                'YLim',[-1.1 0.11],...
                'ZLim',[-0.01 1.01]...
             )
            xlabel('x') 
            ylabel('y')
            zlabel('z')
        end
        
        function bool = hasBodyForce(obj)
           
            bool = true;
        end
        
        function RHS = bodyForce(obj)
            
            force = obj.getBodyForce();
            
            RHS = -force.data.';
        end
        
        function force = getBodyForce(obj)
            
            mesh = obj.mesh();
            N = mesh.numberOfElements();
            
            density = obj.materialLaw().density();
            data = repmat(density*obj.gravity_, [1,N]);
            density_g = vectorField(data);
            force = mesh.volume() * density_g;
        end
    end
end

