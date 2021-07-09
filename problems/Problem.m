classdef (Abstract) Problem < handle
    %Problem Base class for a simulation case.
    %   Subclass this in order to run a case.
    
    properties
        
        % Simulation mesh.
        mesh_;
        
        % Handles smulation time.
        runTime_;
        
        % General settings parameters of the simulation.
        config_;
        
        % A solid solver
        solidModel_;
        
        % Initial position of vol. and face centroids.
        volX_;
        
        % Initial position for face centroids.
        surfaceX_;
        
        % Callback function called after problem object has received the
        % message informing about solidModel evolution/update.
        solidModelhasEvolvedCallback_;
        
        % If matlab is not running GUI
        runningOnConsole_;
    end
    
    methods
        
        function obj = Problem(setup, params)
            %Problem Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin == 1
                params = [];
            end
            
            % Setup time.
            if setup.nTimeSteps == 0
                
                obj.runTime_ = SteadyTime();
            else
                nTimeSteps = setup.nTimeSteps;
                startTime = setup.startTime;
                endTime = setup.endTime;
                deltaT = (endTime - startTime) / nTimeSteps;
                
                if true
                
                deltaT_o = 0;%deltaT;
                obj.runTime_ = Time(startTime, endTime, deltaT, deltaT_o);
                
                else
                
                deltaT_o = deltaT/2.2;
                obj.runTime_ = Time(startTime, endTime, deltaT, deltaT_o);
                obj.runTime_.increment();
                obj.runTime_.resetTimeIndex();
                obj.runTime_.pushDeltaT(deltaT_o);
                obj.runTime_.pushDeltaT(deltaT_o);
                end
            end
            
            obj.mesh_ = ...
                fvMesh...
                (...
                    obj.runTime_,...
                    fvmReadOpenFoamMesh(setup.casePath)...
                );
            
            obj.config_ = configGenerator();
            
            obj.setInitialConfiguration();
            
            obj.solidModelhasEvolvedCallback_ = ...
                tryGetOrDefault...
                (...
                    params,...
                    'solidModelhasEvolvedCallback',...
                    @(the_problem) the_problem...
                );
            
            obj.setNumberOfCorrectors(...
                tryGetOrDefault(params, 'numberOfCorrectors', 100));
            
            obj.setRelaxation(...
                tryGetOrDefault(params, 'relaxation', 1.0));
            
            obj.setSolutionTolerance(...
                tryGetOrDefault(params, 'solutionTolerance', 1e-6));
            
            obj.setAlternativeTolerance(...
                tryGetOrDefault(params, 'alternativeTolerance', 1e-7));
            
            obj.runningOnConsole_ = ...
                tryGetOrDefault(params, 'runningOnConsole', false);
        end
        
        %% BEG - Data member access
        function mesh = mesh(obj)
            % mesh Returns the simulation mesh.
            
            mesh = obj.mesh_;
        end
        
        function obj = setMesh(obj, mesh)
            
           obj.mesh_ = mesh; 
        end
        
        function runTime = runTime(obj)
            % runTime Returns the object which handles smulation time.
            
            runTime = obj.runTime_;
        end
        
        function config = configuration(obj)
            %configuration General settings parameters of the simulation
            
            config = obj.config_;
        end
        
        function solidModel = solidModel(obj)
            % solidModel() return the solid model.
            
            solidModel = obj.solidModel_;
        end
        
        function X = volX(obj)
            %volX Initial position.
            
            X = obj.volX_;
        end
        
        function setVolX(obj, newValue)
            %setVolX Initial position
            %   setVolX(newValue) set the initial position using
            %   the volTensorField newValue
            
            obj.volX_ = newValue;
        end
        
        function X = surfaceX(obj)
            %surfaceX Initial position.
            
            X = obj.surfaceX_;
        end
        
        function setSurfaceX(obj, newValue)
            %setSurfaceX Initial position
            %   setSurfaceX(newValue) set the new initial position using
            %   the surfaceTensorField newValue.
            
            obj.surfaceX_ = newValue;
        end
        %% END
        
        %% BEG - Configuration access
        function nCorrectors = numberOfCorrectors(obj)
           % numberOfCorrectors % Returns the number of momentum 
           % correctors.
            
           nCorrectors = obj.config_.nCorrectors;
        end
        
        function obj = setNumberOfCorrectors(obj, value)
           % numberOfCorrectors % Maximum number of momentum correctors.
            
           obj.config_.nCorrectors = value;
        end
        
        function obj = applyUnderRelaxation(obj, value)
            % applyRelaxation Turn on linear system relaxation.
            
            obj.config_.applyUnderRelaxation = value;
        end
        
        function relaxation = relaxation(obj)
           % relaxation Returns the linear system relaxation.
           
            relaxation = obj.config_.relaxation;
        end
        
        function obj = setRelaxation(obj, value)
            % setRelaxation Sets the linear system relaxation.
            
            obj.config_.relaxation = value;
        end
        
        function obj = setSolutionTolerance(obj, value)
            % setSolutionTolerance Sets the linear momentum convergence 
            % tolerance.
            
            obj.config_.solutionTolerance = value;
        end
        
        function obj = setAlternativeTolerance(obj, value)
            % setAlternativeTolerance Sets the linear momentum convergence 
            % tolerance.
            
            obj.config_.alternativeTolerance = value;
        end
        
        function obj = setLinearSystemTolerance(obj, value)
            % setLinearSystemTolerance Sets the linear system tolerance.
            
            obj.config_.linearSystemTolerance = value;
        end
        
        function obj = setAmgRelaxation(obj, value)
            % setAmgRelaxation Set the AMG relaxation parameter.
            
            obj.config_.amgRelaxation = value;
        end
        %% END
        
        function setInitialConfiguration(obj)
            %setInitialConfiguration Set the initial configuration.
            %   The default is to set the solid as undeformed.
            
            obj.setVolX(obj.mesh_.r_C());
            obj.setSurfaceX(obj.mesh().r_Cf());
        end
        
        function obj = printTime(obj)
            %printTime Print to console a formated time.
            
            fprintf('Time = %d\n', obj.runTime_.value());
        end
        
        function setCamera(obj)
            % Empty
        end
        
        function obj = statistics(obj, solidModel) %#ok<INUSD>
            % Empty
        end
        
        function bool = hasBodyForce(obj)
            % hasBodyForce Override this function if you need to use
            %   body force. Change the sentence below to "bool = true", and
            %   implement the function RHS = bodyForce(obj).
            
            bool = false;
        end
        
        function bool = includeAcceleration(obj)
            % includeAcceleration Override this function if you need to
            %   remove the temporal term acceleration. Change the sentence
            %   below to "bool = false".
            
            bool = true;
        end
        
        function setSolidModel(obj, params)
            % setSolidModel(params) set up a solid model using params.
            %   As default: Rubber, Segregated and Hooke
            
            materialEnum = ...
                tryGetOrDefault...
                (...
                    params,...
                    'material',...
                    materials.List.Rubber...
                );
            
            solidModelSolverEnum = ...
                tryGetOrDefault...
                (...
                    params,...
                    'solidModelSolver',...
                    solidModelSolvers.List.Segregated...
                );
            
            materialLawEnum = ...
                tryGetOrDefault...
                (...
                    params,...
                    'materialLaw',...
                    materialLaws.List.HookeanElastic ...
                );
            
            material = obj.createMaterial(materialEnum);
            materialLaw = obj.createMaterialLaw(materialLawEnum, material);
            
            mechanicalModel = ...
                obj.createMechanicalModel...
                (...
                    materialLaw,...
                    solidModelSolverEnum...
                );
        
            solidModelSolver = ...
                obj.createSolidModelSolver...
                (...
                    solidModelSolverEnum,...
                    mechanicalModel...
                );

             obj.solidModel_ = solidModelSolver;
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
                    varargin{:}...
                );
        end
        
        function mechanicalModel = createMechanicalModel...
                (...
                    obj,...
                    materialLaw,...
                    solidModelSolverEnum,...
                    varargin...
                )
            
            mechanicalModel = ...
                mechanicalModels.createMechanicalModel...
                (...
                    materialLaw,...
                    solidModelSolverEnum,...
                    varargin{:}...
                );
        end
        
        function solidModelSolver = ...
                createSolidModelSolver...
                (...
                    obj,...
                    solidModelSolverEnum,...
                    mechanicalModel,...
                    varargin...
                )
            
            solidModelSolver = ...
                solidModelSolvers.createSolidModelSolver...
                (...
                    solidModelSolverEnum,...
                    mechanicalModel,...
                    obj,...
                    varargin{:}...
                );
        end
        
        function materialLaw = materialLaw(obj)
            
            materialLaw = ...
                obj.solidModel_.mechanicalModel().materialLaw();
        end
        
        function solidModelHasEvolved(obj)
            %solidModelHasEvolved should be called after the solidModel has
            % evolved. Override this method as you wish.
            
            obj.solidModelhasEvolvedCallback_(obj);
        end
        
        function obj = solve(obj)

            solidModel = obj.solidModel();

            %if ~obj.runningOnConsole_

                suffix = strcat('-', string(solidModel.problem().runTime().timeIndex()));
                obj.meshDrawer().run(...
                    solidModel.volU(), solidModel.volGradU());
                obj.setCamera();
                drawnow;
            %end

            printHeader();

            while obj.runTime().run()

                obj.runTime().increment();

                obj.printTime();

                solidModel.evolve();

                obj.solidModelHasEvolved();

                %if ~obj.runningOnConsole_
                    suffix = strcat('-', string(solidModel.problem().runTime().timeIndex()));
                    obj.meshDrawer().run(...
                        solidModel.volU(), solidModel.volGradU());
                %end
            end

            printFooter();
        end
    end
    
    methods (Abstract)
       
        meshDrawer = meshDrawer(obj);
    end
    
    methods (Static)
       
        function namesList = getNameOfCases(classPrefix)

            filePath = Problem.getRootFolder(classPrefix);
            list = dir([filePath filesep 'meshes']);
            namesList = {list(3:end).name};
        end
        
        function filePath = getRootFolder(classPrefix)

            filePath = fileparts(which(classPrefix));

            if isempty(filePath)
                error('Could not find class prefix: %s\n', classPrefix);
            end
        end  
        
        function fullPaths = getMeshPaths(fullPath, className)
            %getMeshPaths return a cell with all mesh paths.
        
            root = fullPath(1:end-length(className));
            root = [root , 'meshes'];
            list = dir(root);
            folders = {list(3:end).name};
            fullPaths = cell(1,length(folders));
            
            for i = 1:length(folders)
                
                fullPaths{i} = [ root, filesep , folders{i} ];
            end
        end
    end
end

