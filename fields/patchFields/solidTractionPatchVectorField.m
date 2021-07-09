classdef solidTractionPatchVectorField < patchVectorField
    %solidTractionPatchVectorField Summary of this class goes here
    %   The solidTraction boundary condition is implemented using 
    %   options 2: 1) and 3) in [Pag. 278 - Moukalled] for internal faces
    %   and using taylor series to calculate the value of phi at the face's
    %   centroid.
    %
    % TODO: 
    %  - To be compatible with solids4Foam:
    % NAME
    % {
    %     type            solidTraction;
    %     traction        uniform (0 0 0);
    %     pressure        uniform 0;
    %     value           uniform (0 0 0);
    % }
    
    properties
        
        % External traction acting over the boundary.
        externalTraction_;
    end

    methods
        
        function obj = ...
                solidTractionPatchVectorField...
                (...
                    mesh,...
                    boundary,...
                    inputVectorField...
                 )
            %solidTractionPatchVectorField Construct an instance of 
            %   this class. Detailed explanation goes here
            
            args{1} = mesh;
            args{2} = boundary.index;
            
            if nargin == 3
                
                args{3} = inputVectorField;
            end
            
            obj@patchVectorField(args{:});
            
            % Does boundary have traction value?
            if isfield(boundary, 'traction')
               
                value = boundary.traction;
                obj.externalTraction_ = patchVectorField(args{1:2}, value);
            else
                % Set traction to zero.
                obj.externalTraction_ = patchVectorField(args{1:2});
            end
        end
        
        function externalTraction = externalTraction(obj)
            
            externalTraction = obj.externalTraction_;
        end
        
        function obj = setExternalTraction(obj, externalTraction)
           
            obj.externalTraction_ = externalTraction;
        end
        
        function obj = update(obj)
            
            % Nothing to do.
        end
        
        function obj = updateDivImpKgradUfluxes(obj, solidModel)
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            divImpKgradU = ...
                solidModel.tractionBoundaryLaplacianOfimpKGradU...
                (...
                    obj.externalTraction_,...
                    obj...
                );
            
            divImpKgradU_LHS_ = [];
            divImpKgradU_RHS_ = -(divImpKgradU.field().data).';
            
            obj = ...
                obj.setDivImpKgradUfluxes...
                ( ...
                    divImpKgradU_LHS_,...
                    divImpKgradU_RHS_...
                );
        end
    end
end

