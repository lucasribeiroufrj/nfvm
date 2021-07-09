classdef Uniaxial < MMSProblem & handle
    %Uniaxial Uniaxial compression/elongation
    %   Case described in [5, Kamojjala - 2013].
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
    %
    % Ref.:
    %   [Kamojjala - 2013] - Verification tests in solid mechanics
    
    properties
        
        % Compression factor
        Lambda_;
        
        % The time derivative of the deformation tensor.
        dFdt_;
        
        % To draw to mesh
        meshDrawer_;
        
        % Boundary conditions for displacement U 
        boundaryConditions_;
        
        % Bool to indicate whether to use only Dirichlet boundary cond. or
        % not.
        useTraction_;
        
        % Bool to indicate whether to use symmetry boundary cond. or
        % not (perpendicular-to-compressed-face boundaries).
        useSymmetryBoundaries_;
    end
    
    methods
        
        function obj = Uniaxial(casePath, params)
            
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
            
            obj.Lambda_ = tryGetOrDefault(params, 'Lambda', 0.5);
            
            obj.setNumberOfCorrectors(...
                tryGetOrDefault(params, 'numberOfCorrectors', 50));
            
            obj.useTraction_ = ...
                tryGetOrDefault(params, 'useTraction', false);
            
            obj.useSymmetryBoundaries_ = ...
                tryGetOrDefault...
                (...
                    params,...
                    'useSymmetryBoundaries',...
                    true...
                );
            
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
            
            mesh = obj.mesh().spaceMesh();
            
            nTotalElements = mesh.numberOfElements...
                + mesh.numberOfBElements;
            M = zeros(3,3);
            M(1,1) = obj.Lambda_ - 1;
            dFdt = repmat(M, [1,1,nTotalElements]);
            obj.dFdt_ = volTensorField(mesh, tensorField(dFdt));
        end
        
        function obj = setBoundaryConditions(obj)
            
            mesh = obj.mesh();
            
            for iBoundary = 1:mesh.numberOfBoundaries()
                
                boundary = mesh.boundaries(iBoundary);
                bName = boundary.userName;
                
                if strcmp(bName,'left')

                    kind = 'MMSDisplacement';
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
        
        function value = deltaLength(obj, time, iBoundary)
            
            mesh = obj.mesh();
            nf = mesh.nf().boundaryField(3);
            solidModel = obj.solidModel();
            
            value = obj.analyticalVolU(time).boundaryField(iBoundary).field;
            return;
            
            % Is incremental approach?
            if false && strfind(lower(class(solidModel)), ...
                    lower('BlockCoupledWithNewBoundary'))
                
                if mesh.timeMesh.isFirstTimeStep()
                    t_o = mesh.timeMesh.startTime();
                else
                    t_o = mesh.timeMesh.valueOld();
                end
                t_ = mesh.timeMesh.value();

                x   = 1 + (obj.Lambda() - 1)*t_;
                x_o = 1 + (obj.Lambda() - 1)*t_o;
                du = x - x_o;
            
                value = vectorField(du*nf);
            else
                t_ = mesh.timeMesh.value();
                x = (1 + (obj.Lambda() - 1)*t_);
                x_o = 1;
                du = x - x_o;
                
                value = vectorField(du*nf);
            end
        end
        
        function obj = setMeshDrawer(obj)
           
            obj.meshDrawer_ = BoundaryMeshDrawer(obj.mesh());
        end
        
        function U = U(obj)
            
           U = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function U_o = U_o(obj)

            deltaT_o = obj.runTime.deltaT_o();
            V_oo = obj.dFdt_*obj.volX();
            U_o = obj.U_oo() + V_oo*deltaT_o;
        end
        
        function U_oo = U_oo(obj)
            
            % Displacement is set to zero for now.
            U_oo = volVectorField(obj.mesh());
        end
        
        function DU = DU(obj)
           
            DU = obj.U() - obj.U_o();
        end
        
        function Lambda = Lambda(obj)
           
            Lambda = obj.Lambda_;
        end
        
        %function F = computeTensorFieldF(obj, time, numberOfElements)
        function F = computeTensorFieldF(obj, time, positions)
            %computeTensorFieldF Computation of the deformation F.
            %   computeTensorFieldF(time, numberOfElements) returns F as
            %   tensorField at a specific time.
            nElements = positions.numberOfElements;
            M = eye(3,3);
            M(1,1) = 1 + (obj.Lambda_ - 1)*time;
            data = repmat(M, [1,1,nElements]);
            F = tensorField(data);
        end
        
        function meshDrawer = meshDrawer(obj)
            
            meshDrawer = obj.meshDrawer_;
        end
        
        function setCamera(obj)
            
            dx = 0;
            if obj.Lambda_> 1
                dx = obj.Lambda_ - 1;
            end
            
            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],...
                'XLim',[-0.1 1.1 + dx],...
                'YLim',[-0.1 1.1],...
                'ZLim',[-0.1 1.1]...
             )
            xlabel('x') 
            ylabel('y')
            zlabel('z')
            
            camorbit(20, 20,'data',[0 1 0]);
        end
        
        function x = computeVectorField_x(obj, time, X)
            %computeVectorField_x Computation of the position x.
            %   computeVectorField_x(time, iX) returns x as
            %   vectorField at a specific time. 
            %   position: must be a vectorField.
            
            F = obj.computeTensorFieldF(time, X);
            x = F*X;
        end
        
        function bool = includeAcceleration(obj)
            % includeAcceleration Override this function if you need to
            %   remove the temporal term acceleration. Change the sentence
            %   below to "bool = false".
            
            bool = false;
        end
    end
end

