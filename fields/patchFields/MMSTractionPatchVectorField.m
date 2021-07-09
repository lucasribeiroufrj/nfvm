classdef MMSTractionPatchVectorField < solidTractionPatchVectorField
    %MMSTractionPatchVectorField Boundary condition for MMS* test cases.
    %   Detailed explanation goes here
    %
    % MMS - Method of Manuctured solutions.
    
    properties
        
        % Problem reference
        problem_
        
        % timeIndex
        timeIndex_;
        
        % A material law
        materialLaw_;
    end

    methods
        
        function obj = ...
                MMSTractionPatchVectorField...
                (...
                    mesh,...
                    boundary,...
                    inputVectorField...
                )
            %MMSTractionPatchVectorField Construct an instance of 
            %   this class. Detailed explanation goes here

            args{1} = mesh;
            args{2} = boundary;
            
            if nargin == 3
                
                args{3} = inputVectorField;
            end
                
            obj@solidTractionPatchVectorField(args{:});
            
            obj.problem_ = boundary.problem;
            obj = obj.initialize();
        end
        
        function obj = initialize(obj)
            
            obj.timeIndex_ = -1;
        end
        
        function obj = update(obj)
           
            timeMesh = obj.mesh().timeMesh();
            currentTimeIndex = timeMesh.timeIndex();
            
            % A new time step?
            if obj.timeIndex_ ~= currentTimeIndex
                
                obj.timeIndex_ = currentTimeIndex;
            else % Already updated the boundary condition?
                return;
            end
            
            idx = obj.index();
            nf = obj.mesh().nf().boundaryField(idx);
            solidModel = obj.problem_.solidModel();
            mechanicalModel = solidModel.mechanicalModel();
            
            time = obj.problem_.runTime().value();
            gradU = obj.problem_.analyticalVolGradU(time)...
                .boundaryField(idx);
            P = mechanicalModel.computePiolaField(gradU);
            externalTraction = P*nf;

            obj = obj.setExternalTraction(externalTraction);
        end

        function obj = updateDivImpKgradUfluxes(obj, solidModel)
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            obj = obj.update();
            
            obj = ...
                obj.updateDivImpKgradUfluxes@solidTractionPatchVectorField...
                (...
                    solidModel...
                );
        end
    end
end

