classdef emptyPatchVectorField < patchVectorField
    %zeroDisplacementPatchVectorField Summary of this class goes here
    %   Detailed explanation goes here
    %
    % TODO: 
    %  - To be compatible with solids4Foam:
    % NAME
    % {
    %     type            fixedDisplacement;
    %     displacementSeries
    %     {
    %         fileName        "$FOAM_CASE/constant/timeVsDisp";
    %         outOfBounds     clamp;
    %     }
    %     value           uniform (0 0 0);
    % }
    
    properties
        
    end
    
    methods
        
        function obj = ...
                emptyPatchVectorField...
                (...
                    mesh,...
                    boundary...
                 )  % inputVectorField...
                %)
            %emptyPatchVectorField Construct an instance of 
            %
            %   emptyPatchVectorField(mesh, boundary) is a patch
            %   with geometric boundary information given by 
            %   mesh.bondaries(boundary.index).
            %
            %   emptyPatchVectorField(mesh, index, 
            %   inputVectorField) is the same as above plus field values.
            %
            
            args{1} = mesh;
            args{2} = boundary.index;
            
            if nargin == 3
                
                args{3} = inputVectorField;
            end
            
            obj@patchVectorField(args{:});
        end
        
        function obj = updateDivImpKgradUfluxes(obj, solidModel) %#ok<INUSD>
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            % Nothing to do.
        end
        
        function obj = correct(obj, volU, volGradU, surfaceGradU) %#ok<INUSD>
            
            % Nothing to do.
        end
    end
end

