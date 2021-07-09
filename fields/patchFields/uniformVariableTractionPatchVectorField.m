classdef uniformVariableTractionPatchVectorField ...
        < solidTractionPatchVectorField
    %uniformVariableTractionPatchVectorField Traction boundary condition.
    %   Uniform means same traction value for all faces, but this value
    % may be a function of time. This traction is along the face's
    % normal.
    
    properties
        
        % A handle to a vector-valued function traction=traction(time)
        tractionFunction_;
        
        % timeIndex
        timeIndex_;
        
        % Problem reference
        problem_;
    end
    
    methods
        function obj = ...
                uniformVariableTractionPatchVectorField...
                (...
                    mesh,...
                    boundary,...
                    inputVectorField...
                )
            
            args{1} = mesh;
            args{2} = boundary;
            
            if nargin == 3
                
                args{3} = inputVectorField;
            end
            
            obj@solidTractionPatchVectorField(args{:});
            
            obj.tractionFunction_ = boundary.tractionFunction;
            obj.problem_ = boundary.problem;
            
            obj = obj.initialize();
        end
        
        function obj = initialize(obj)
            
            obj.timeIndex_ = -1;
        end
        
        function obj = updateDivImpKgradUfluxes(obj, solidModel)
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            obj = obj.correct();
            obj = ...
                obj.updateDivImpKgradUfluxes@solidTractionPatchVectorField...
                (...
                    solidModel...
                );
        end
        
        function obj = correct(obj, volU, volGradU, surfaceGradU) %#ok<INUSD>
            %correct Computes the Neumann (traction) boundary at time t.
            %   It will be updated only once for each time step.
            
            timeMesh = obj.mesh().timeMesh();
            currentTimeIndex = timeMesh.timeIndex();
            
            if obj.timeIndex_ ~= currentTimeIndex
                
                obj.timeIndex_ = currentTimeIndex;
            else
                return;
            end
            
            currentTime = timeMesh.value();
            traction = obj.tractionFunction_(currentTime);
            
            obj = obj.setExternalTraction(traction);
        end
    end
    
end

