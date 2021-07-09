classdef calculatedPatchScalarField < patchScalarField
    %patchScalarField An abstraction of a boundary condition.
    %   Detailed explanation goes here
    %
    % TODO:
    % - We should use only a sub-mesh, but should also have a ref. to the
    % main mesh.
    
    properties

    end
    
    methods
        
        function obj = calculatedPatchScalarField(mesh, index)
            %patchScalarField Construct an instance of this class
            %   Detailed explanation goes here
            
            obj@patchScalarField(mesh, index)
        end
        
        function output = field(obj,i,j)
            
            if nargin == 1
            
                output = obj;
                
            elseif nargin == 2

                output = obj.subsref(i);

            elseif nargin == 3

                output = obj.subsref(i,j);

            else
                error('Use one or two parameters');
            end
        end
    end
end

