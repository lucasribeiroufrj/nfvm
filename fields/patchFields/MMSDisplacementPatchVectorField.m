classdef MMSDisplacementPatchVectorField < patchVectorField
    %MMSDisplacementPatchVectorField Boundary condition for MMS* test cases.
    %   Detailed explanation goes here
    %
    % MMS - Method of Manuctured solutions.
    
    properties
        
        % Problem reference
        problem_
        
        % timeIndex
        timeIndex_;
    end
    
    methods
        function obj = ...
                MMSDisplacementPatchVectorField...
                (...
                    mesh,...
                    boundary...
                )
            %MMSDisplacementPatchVectorField 
            %   Detailed explanation goes here
            
            args{1} = mesh;
            args{2} = boundary.index;
            
            obj@patchVectorField(args{:});
            
            obj.problem_ = boundary.problem;
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
            
            % A new time step?
            if obj.timeIndex_ ~= currentTimeIndex
                
                obj.timeIndex_ = currentTimeIndex;
            else % Already updated the boundary condition?
                return;
            end
            
            iBoundary = obj.index();
            time = timeMesh.value();
            U = obj.problem_.analyticalVolU(time).boundaryField(iBoundary);
            obj = obj.setVectorField(U.field());
        end
        
        function obj = correct(obj, volU, volGradU, surfaceGradU) %#ok<INUSD>
            
            obj = obj.update();
        end
    end
end
