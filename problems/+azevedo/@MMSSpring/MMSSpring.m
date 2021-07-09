classdef MMSSpring < MMSProblem & handle
    %MMSSpring Case similar to [imechanica] with analytical answer 
    %   (Method of Manufactured Solutions). 
    %
    % TODO:
    %   - 
    %
    % Nomenclature.: 
    % - The letters U stand for displacement vector
    % respectively.
    % - The U_o, U_oo stands for U at old, and old old time respectively.
    %
    % Ref.:
    % [imechanica] - https://imechanica.org/node/1357
    
    properties
        
        % To draw to mesh
        meshDrawer_;
        
        % Boundary conditions for displacement U 
        boundaryConditions_;
        
        % Phase
        w_;
    end
    
    methods
        
        function obj = MMSSpring(casePath, params)
            
            nTimeSteps = ...
                tryGetOrDefault(params, 'numberOfTimeSteps', 5);
            
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                    'nTimeSteps', nTimeSteps,...
                    'startTime', 0.0,...
                    'endTime', 1.0...
                );
            
            obj@MMSProblem(setup, params);
            
            obj.w_ = tryGetOrDefault(params, 'phase', pi/2);
            
            obj.setFields();
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
        
        function obj = setFields(obj)
            
            %mesh = obj.mesh().spaceMesh();
        end
        
        function obj = setBoundaryConditions(obj)
            
            for iBoundary=1:obj.mesh().numberOfBoundaries
                
                boundary = obj.mesh().boundaries(iBoundary);
                bName = boundary.userName;
                boundCondition = [];
                
                if strcmp(bName,'front') || strcmp(bName,'back') ...
                        || strcmp(bName,'left') || strcmp(bName,'right') ...
                        || strcmp(bName,'top') || strcmp(bName,'bottom')

                    boundCondition.kind = 'MMSDisplacement';
                    boundCondition.problem = obj;
                else
                    ME = MException('MMSSpring:boundaryCondition',...
                        'Invalid boundary name: %s', bName);
                    throw(ME)
                end
                
                boundCondition.index = iBoundary;
                obj.boundaryConditions_{iBoundary} = boundCondition;
            end
        end
        
        function bool = hasBodyForce(obj)
           
            bool = true;
        end
        
        function RHS = bodyForce(obj)
            
            force = obj.getBodyForce();
            
            RHS = -force.data.';
        end
        
        function force = getBodyForce(obj)
            
            if ~isa(obj.materialLaw(), 'materialLaws.NeoHookean')
            
                error('Only NeoHookean is supported yet');
            end
            
            mesh = obj.mesh();
            
            X = obj.volX().internalField();
            nElements = X.numberOfElements;
            data = zeros(3,nElements);
            dataX = X.data;
            density = obj.materialLaw().density();
            mu = obj.materialLaw().material().mu();
            lambda = obj.materialLaw().material().lambda();
            w = obj.w_;
            
            for iElement = 1:nElements
                
                X1 = dataX(1,iElement);
                t = mesh.timeMesh.value();
                dF11dX1 = 2*sin(w*t);
                d2u1dX12 = 2*sin(w*t);
                du1dX1 = X1*d2u1dX12;
                dinvFT11dX1 = -(1 + du1dX1)^(-2)*(d2u1dX12);
                J = 1 + 2*X1*sin(w*t);
                invFT11 = 1.0/(1 + du1dX1);
                
                d2u1dt2 = -(X1*w)^2*sin(w*t);
                f1 = mu*(dF11dX1 - dinvFT11dX1) ...
                    + lambda*(1.0/J*invFT11 + log(J)*dinvFT11dX1);
                
                data(:,iElement) = [d2u1dt2 - f1/density; 0; 0];
            end
            
            force = vectorField(data)*mesh.volume()*density;
        end
        
        function obj = setMeshDrawer(obj)
           
            obj.meshDrawer_ = BoundaryMeshDrawer(obj.mesh());
        end
        
        function DU = DU(obj)
           
            DU = obj.U() - obj.U_o();
            error('Not implemented!');
        end
        
        function F = computeTensorFieldF(obj, time, numberOfElements)
            %computeTensorFieldF Computation of the deformation F.
            %   computeTensorFieldF(time, numberOfElements) returns F as
            %   tensorField at a specific time.
            
            error('Not implemented!');
            M = eye(3,3);
            M(1,2) = obj.shearFactor_*time;
            F = tensorField(numberOfElements);
            F(1:end) = M;
        end
        
        function meshDrawer = meshDrawer(obj)
            
            meshDrawer = obj.meshDrawer_;
        end
        
        function setCamera(obj)
            
            xLim = 2.2;
            
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
        
        function x = computeVectorField_x(obj, time, X)
            %computeVectorField_x Computation of the position x.
            %   computeVectorField_x(time, iX) returns x as
            %   vectorField at a specific time. 
            %   position: must be a vectorField.
            
            nElements = X.numberOfElements;
            data = zeros(3,nElements);
            dataX = X.data;
            w = obj.w_;
            
            for iElement = 1:nElements
                
                X1 = dataX(1,iElement);
                X2 = dataX(2,iElement);
                X3 = dataX(3,iElement);
                U1 = X1*X1*sin(w*time);
                x = [...
                        X1 + U1;
                        X2;
                        X3...
                    ];
                
                data(:,iElement) = x;
            end
            
            x = vectorField(data);
        end
    end
end

