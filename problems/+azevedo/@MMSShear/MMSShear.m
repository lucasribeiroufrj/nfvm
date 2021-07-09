classdef MMSShear < MMSProblem & handle
    %MMSShear Shear case with analytical answer (Method of Manufactured
    %   Solutions).
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
    % [Bonet - 2008] - Nonlinear continuum mechanics for finite element
    % analysis
    
    properties
        
        % x_1 = X_1 + shearFactor_*X_2.
        % x_2 = X_2.
        % x_3 = X_3.
        shearFactor_;
        
        % The time derivative of the deformation tensor.
        dFdt_;
        
        % To draw to mesh
        meshDrawer_;
        
        % Boundary conditions for displacement U 
        boundaryConditions_;
        
        % A vol vector field V | V(i) = [X2; 0, 0], X2 being the 2nd
        % coordinate of FV's centroid i.
        volShearDirection_;
        
        % Bool to indicate whether to use only Dirichlet boundary cond. or
        % not.
        useTraction_
    end
    
    methods
        
        function obj = MMSShear(casePath, params)
            
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
            
            obj@MMSProblem(setup, params);
            
            obj.shearFactor_ = ...
                tryGetOrDefault(params, 'shearFactor', 0.45);
            
            obj.setNumberOfCorrectors(...
                tryGetOrDefault(params, 'numberOfCorrectors', 50));
            
            obj.useTraction_ = ...
                tryGetOrDefault(params, 'useTraction', false);
            
            obj.setFields();
            obj.setBoundaryConditions();
            obj.setMeshDrawer();
            
            obj.setSolidModel(params);
            
            obj.computeShearDirection();
        end
        
        function obj = setFields(obj)
            
            mesh = obj.mesh().spaceMesh();
            
            nTotalElements = mesh.numberOfElements...
                + mesh.numberOfBElements;
            M = zeros(3,3);
            M(1,2) = obj.shearFactor_;
            dFdt = repmat(M, [1,1,nTotalElements]);
            obj.dFdt_ = volTensorField(mesh, tensorField(dFdt));
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
                
                if strcmp(bName,'bottom')

                    %kind = 'MMSDisplacement';
                    kind = 'fixedDisplacement';
                else
                    if obj.useTraction_

                        kind = 'MMSTraction';
                    else
                        kind = 'MMSDisplacement';
                    end
                end
                
                obj.boundaryConditions_{iBoundary}.index = iBoundary;
                obj.boundaryConditions_{iBoundary}.problem = obj;
                obj.boundaryConditions_{iBoundary}.kind = kind;
            end
        end
        
        function value = displacementAtBoundary(obj, time, iBoundary)

            mesh = obj.mesh();
            solidModel = obj.solidModel();
            t_ = time;
            shearDir = obj.volShearDirection_.boundaryField(iBoundary);
            
            % Is incremental approach?
            if false && strfind(lower(class(solidModel)), ...
                    lower('BlockCoupledWithNewBoundary'))
                
                if mesh.timeMesh.isFirstTimeStep()
                    
                    t_o = mesh.timeMesh.startTime();
                else
                    t_o = mesh.timeMesh.valueOld();
                end
                
                value = (t_ - t_o)*obj.shearFactor_*shearDir;
            else
                value = t_*obj.shearFactor_*shearDir;
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
            
            deltaT_o = obj.runTime().deltaT_o;
            V_oo = obj.dFdt_*obj.volX();
            U_o = obj.U_oo() + V_oo*deltaT_o;
        end
        
        function U_oo = U_oo(obj)
            %Uoo Displacement field at time "old-old"
            
            U_oo = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function DU = DU(obj)
           
            DU = obj.U() - obj.U_o();
        end
        
        function F = computeTensorFieldF(obj, time, X)
            %computeTensorFieldF Computation of the deformation F.
            %   computeTensorFieldF(time, X) returns F as
            %   tensorField at specific time and positions given by the
            %   vectorField X.
            
            M = eye(3,3);
            M(1,2) = obj.shearFactor_*time;
            data = repmat(M, [1,1,X.numberOfElements]);
            F = tensorField(data);
        end
        
        function x = computeVectorField_x(obj, time, X)
            %computeVectorField_x Computation of the position x.
            %   computeVectorField_x(time, iX) returns x as
            %   vectorField at specific time and 
            %   position: must be a vectorField.
            
            F = obj.computeTensorFieldF(time, X);
            x = F*X;
        end
        
        function meshDrawer = meshDrawer(obj)
            
            meshDrawer = obj.meshDrawer_;
        end
        
        function setCamera(obj)
            
            xLim = 1.2;
            
            if obj.shearFactor_ > 0
                xLim = 1 + obj.shearFactor_ + 0.1;
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
        
        function computeShearDirection(obj)
            %computeShearDirection set volShearDirection_. A volVectorField
            %   V | V(i) = [X2; 0, 0], X2 being the 2nd
            %   coordinate of FV's centroid i.
            %   Override it if you want another direction.
            
            volX = obj.volX();
            e1 = obj.mesh.volE1();
            e2 = obj.mesh.volE2();
            obj.volShearDirection_ = (volX^e2)*e1;
        end
        
        function bool = includeAcceleration(obj)
            % includeAcceleration Override this function if you need to
            %   remove the temporal term acceleration. Change the sentence
            %   below to "bool = false".
            
            bool = false;
        end
    end
end

