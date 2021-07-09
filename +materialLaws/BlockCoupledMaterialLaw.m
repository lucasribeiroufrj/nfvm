classdef (Abstract) BlockCoupledMaterialLaw ...
        < materialLaws.MaterialLaw & handle
    % MaterialLaw Interface class for a BC Material law.
    %   Subclass this if you want implement a material law that conforms to
    %   the BlockCoupled solver.
    
    methods (Abstract)
        
        % mu Shear modulus (also known as G).
        mu = mu(obj)
        
        % lambda Lame's first parameter.
        lambda = lambda(obj)    
    end
end

