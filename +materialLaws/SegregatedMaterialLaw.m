classdef (Abstract) SegregatedMaterialLaw ...
        < materialLaws.MaterialLaw & handle
    % MaterialLaw Interface class for a BC Material law.
    %   Subclass this if you want implement a material law that conforms to
    %   the Segregated solver.
    
    methods (Abstract)
                
        % Computes the first-piola stress tensor for a given displacement
        % gradient field. The gradU parameter can be a volume or a surface
        % fields.
        Piola = computePiolaField(obj, gradU)
        
        % Returns the implicit stiffness field as a surfaceField. The
        % function computePiolaField should have been called for an
        % updated impK.
        surfaceImpK = surfaceImpK(obj)
        
        % Returns the implicit stiffness field as a volField. The
        % function computePiolaField should have been called for an
        % updated impK.
        volImpK = volImpK(obj)
    end
end

