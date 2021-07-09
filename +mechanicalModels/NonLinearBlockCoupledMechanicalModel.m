classdef NonLinearBlockCoupledMechanicalModel...
        < mechanicalModels.MechanicalModel & handle
    %NonLinearBlockCoupledMechanicalModel Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        % The MaterialLaw. Might be HookeanElastic, OgdenStoraker, etc.
        materialLaw_;
    end
    
    methods
        function obj = NonLinearBlockCoupledMechanicalModel(materialLaw)
            %NonLinearBlockCoupledMechanicalModel Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.materialLaw_ = materialLaw;
        end
        
        % Beg - mechanicalModels.MechanicalModel interface implementation
        %
        function Piola = computePiolaField(obj, gradU)
            % computePiolaField computes the first-piola stress tensor for
            %   a given displacement gradient field. The gradU parameter
            %   can be a volume or a surface fields.
            
            Piola = obj.materialLaw().computePiolaField(gradU);
        end
        
        function materialLaw = materialLaw(obj)
            % materialLaw Returns the material law
           
            materialLaw = obj.materialLaw_;
        end
        
        function density = density(obj)
            % density returns the density
            
            density = obj.materialLaw().density();
        end
        % End - mechanicalModels.MechanicalModel interface implementation
        
        function Td = computeSurfaceTd(obj, F, d)
            % computeSurfaceTd return M . n (the . at 2). T
            % The deformation gradient parameter F must be a surface
            % field. The parameter d is the row used from d.
            
            Td = obj.materialLaw().computeSurfaceTd(F,d);
        end

        function Sigma = computeSigmaField(obj, gradU)
            % computeSigmaField computes the second-piola stress tensor for
            %   a given displacement gradient field. The gradU parameter
            %   can be a volume or a surface fields.
           
            Sigma = obj.materialLaw().computeSigmaField(gradU);
        end
    end
end

