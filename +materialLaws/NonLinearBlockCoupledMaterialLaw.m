classdef (Abstract) NonLinearBlockCoupledMaterialLaw ...
        < materialLaws.MaterialLaw & handle
    % NLBCMaterialLaw Interface class for NLBC Material law.
    %   Subclass this if you want implement a material law that conforms to
    %   the NonLinearBlockCoupled solver.
    
    methods (Abstract)
        
        % Computes the second-piola stress tensor for a given displacement
        % gradient field. The gradU parameter can be a volume or a surface
        % fields.
        Sigma = computeSigmaField(obj, gradU)
        
        % computeSurfaceTd return M . n (the . at 2).
        % The deformation gradient parameter F must be a surface
        % field. The parameter k is the row used from F.
        Td = computeSurfaceTd(obj, F, k)
    end
end

