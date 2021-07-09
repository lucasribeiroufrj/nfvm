classdef uniformVariableDisplacementPatchVectorField < patchVectorField
    %rigidSphereContactPatchVectorField Normal displament boundary cond.
    %   Uniform means same displacement value for all faces. but this value
    % may be a function of time. This displacemnt is along the face's
    % normal.
    
    properties
        
        % A handle to a scalar-valued function disp=disp(time)
        displacementFunction_;
        
        % timeIndex
        timeIndex_;
    end
    
    methods
        function obj = ...
                uniformVariableDisplacementPatchVectorField...
                (...
                    mesh,...
                    boundary...
                )
            %uniformVariableDisplacementPatchVectorField 
            %   Detailed explanation goes here
            
            args{1} = mesh;
            args{2} = boundary.index;
            
            obj@patchVectorField(args{:});
            
            obj.displacementFunction_ = boundary.displacementFunction;
            
            obj = obj.initialize();
        end
        
        function obj = initialize(obj)
            
            obj.timeIndex_ = -1;
        end
        
        function obj = updateDivImpKgradUfluxes(obj, solidModel)
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            obj = obj.update();
            obj = obj.updateDivImpKgradUfluxes@patchVectorField(solidModel);
        end
        
        function obj = update(obj)
            %correct Computes the dirichlet (displacement) boundary 
            
            timeMesh = obj.mesh().timeMesh();
            currentTimeIndex = timeMesh.timeIndex();
            
            if obj.timeIndex_ ~= currentTimeIndex
                
                obj.timeIndex_ = currentTimeIndex;
            else
                return;
            end
            
            currentTime = timeMesh.value();
            displacement = obj.displacementFunction_(currentTime);
            obj = obj.setVectorField(displacement.field());
        end
    end
end
