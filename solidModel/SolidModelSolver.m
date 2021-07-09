classdef (Abstract) SolidModelSolver < SolidModel
    %SolidModelSolver Derive a classs from this function in order to create
    %   a solid solver.
    
    properties
        
        % Volume displacement field
        volU_;
        
        % TEMPORARY: we should incoporate them into the volU_
        % ATTENTION: Substitute the postfix o_ for oldTime and o__ for
        % oldOldTime.
        volU_o_;
        volU_oo_;
        
        % Gradient of displacement volume field
        volGradU_;
        
        % Gradient of displacement surface field
        surfaceGradU_;
        
        % Handles a material law.
        mechanicalModel_;
        
        % The simulation case.
        problem_;
        
        % Number of empty faces.
        nEmptyFaces_;
    end
    
    methods
        function obj = SolidModelSolver(mechanicalModel, problem)
            % SolidModelSolver Construct an instance of this class
            %   Detailed explanation goes here
            
            obj@SolidModel(problem.mesh());
            
            obj.volU_ = problem.U();
            obj.volU_o_ = problem.U_o();
            obj.volU_oo_ = problem.U_oo();
            obj.volGradU_ = volTensorField(obj.fvMesh());
            obj.surfaceGradU_ = surfaceTensorField(obj.fvMesh());
            obj.mechanicalModel_ = mechanicalModel;
            obj.problem_ = problem;
            obj.setNumberOfEmptyFaces();
        end
        
        %% BEG - Data member access
        function volU = volU(obj)
            
            volU = obj.volU_;
        end
        
        function volU_o = volU_o(obj)
           
            volU_o = obj.volU_o_;
        end
        
        function volU_oo = volU_oo(obj)
           
            volU_oo = obj.volU_oo_;
        end
        
        function volGradU = volGradU(obj)
           
            volGradU = obj.volGradU_;
        end
        
        function surfaceGradU = surfaceGradU(obj)
           
            surfaceGradU = obj.surfaceGradU_;
        end
        
        function mechanicalModel = mechanicalModel(obj)
           
            mechanicalModel = obj.mechanicalModel_;
        end
        
        function problem = problem(obj)
            
            problem = obj.problem_;
        end
        
        function obj = setVolU(obj, newValue)
           
            obj.volU_ = newValue;
        end
        
        function obj = setVolU_o(obj, newValue)
           
            obj.volU_o_ = newValue;
        end
        
        function obj = setVolU_oo(obj, newValue)
           
            obj.volU_oo_ = newValue;
        end
        
        function obj = setVolGradU(obj, newValue)
           
            obj.volGradU_ = newValue;
        end
        
        function obj = setSurfaceGradU(obj, newValue)
           
            obj.surfaceGradU_ = newValue;
        end
        
        function config = config(obj)
            
            config = obj.problem.configuration();
        end
        
        function nEmptyFaces = nEmptyFaces(obj)
            
            nEmptyFaces = obj.nEmptyFaces_;
        end
        %% END
        
        %% BEG - General member functions.
        function obj = setInternalDataVolU(obj, newValue)
            
            obj.volU_ = obj.volU_.setInternalData(newValue);
        end

        function residual = analyticResidual(obj, fieldName)
            %analyticResidual
            %   analyticResidual(fieldName) returns the norm of the
            %   difference between the numerical and the analytical field
            %   fieldName.
            %   fieldName: a string like volU, surfaceU, volGradU, ...
            
            t = obj.fvMesh.timeMesh.value();
            camelName = [upper(fieldName(1)) fieldName(2:end)];
            aGradU = obj.problem_.(['analytical' camelName])(t);
            nGradU = obj.(fieldName)();
            
            % 2-norm
            ignoreBoundary = true;
            residual = fieldNorm(nGradU - aGradU, 2, ignoreBoundary);
            sError = sprintf('%0.2g', residual);
            fprintf(['Analytical residual ' fieldName ': %s\n'], sError);
        end
        
        function obj = setNumberOfEmptyFaces(obj)
           
            boundaryField = obj.volU().boundaryField();
            nEmptyFaces = 0;
            
            for patch = 1:boundaryField.nPatches_
                
                patchVectorField = boundaryField.patchField(patch);
                
                if SolidModelSolver...
                        .isEmptyBoundaryCondition(patchVectorField)
                    
                    nEmptyFaces = nEmptyFaces ...
                        + patchVectorField.numberOfBFaces_;
                end
            end
            
            obj.nEmptyFaces_ = nEmptyFaces;
        end
        
        function obj = correctVolGradU(obj, isToCorrectBoundary)
            % correctVolGradU Caculate volGradU using volU.
            
            if nargin == 1
                isToCorrectBoundary = true;
            end
            
            if isToCorrectBoundary
            
                obj.setVolGradU(...
                    fvc_grad...
                    (...
                        obj.volU(),...
                        obj.fvMesh().Sf(),...
                        obj.config(),...
                        obj.volGradU(),...
                        obj.surfaceGradU()...
                    ));
            else
                obj.setVolGradU(...
                    fvc_grad...
                    (...
                        obj.volU(),...
                        obj.fvMesh().Sf(),...
                        obj.config()...
                    ));
            end
        end
        
        function obj = correctSurfaceGradU(obj)
            % correctSurfaceGradU Interpolates volGradU to surfaceGradU.
            
            % ACHTUNG!!! NO CORRECTION FOR NON-CONJUNCTIONAL MESH.
            obj.setSurfaceGradU(...
                interpolateVolToFaceGrad(obj.volGradU(), obj.volU()));
        end
        %% END - General member functions.
    end
    
    methods(Static)
       
        function bool = boundaryTest(boundaryU, type)
            
            bool = false;
            name = class(boundaryU);
            
            if strfind(lower(name), type)
                
                bool = true;
            end
        end
        
        function bool = isTractionBoundaryCondition(boundaryU)
           
            bool = SolidModelSolver.boundaryTest(boundaryU, 'traction');
        end
        
        function bool = isDisplacementBoundaryCondition(boundaryU)
            
            bool = SolidModelSolver.boundaryTest(boundaryU,'displacement');
        end
        
        function bool = isEmptyBoundaryCondition(boundaryU)
            
            bool = SolidModelSolver.boundaryTest(boundaryU, 'empty');
        end
        
        function bool = isSymmetryBoundaryCondition(boundaryU)
            
            bool = SolidModelSolver.boundaryTest(boundaryU, 'symmetry');
        end
    end
end

