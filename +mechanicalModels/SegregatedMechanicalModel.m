classdef SegregatedMechanicalModel ...
        < mechanicalModels.MechanicalModel & handle
    %SegregatedMechanicalModel Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)

        % The MaterialLaw. Might be Hookean, OgdenStorakerean, etc.
        materialLaw_;
    end

    methods
        function obj = SegregatedMechanicalModel(materialLaw)
            % MechanicalModel Construct an instance of this class
            %   % materialLaw : might be Hookean, OgdenStorakerean, etc.
            
            obj.setMaterialLaw(materialLaw);
        end
        
        function materialLaw = materialLaw(obj)
            % materialLaw Returns the material law
           
            materialLaw = obj.materialLaw_;
        end
        
        function obj = setMaterialLaw(obj, materialLaw)
            % setMaterialLaw sets a material law given a new materialLaw as
            %   input.
            
            obj.materialLaw_ = materialLaw;
        end
        
        function density = density(obj)
            % density returns the density
            
            density = obj.materialLaw().density();
        end
        
        function Piola = computePiolaField(obj, gradU)
            % computePiolaField the first-piola stress tensor for a given
            %   displacement gradient field.
            %   The gradU parameter can be a volume or a surface 
            %   fields.
            
            Piola = obj.materialLaw().computePiolaField(gradU);
        end
        
        function surfaceImpK = surfaceImpK(obj)
            % surfaceImpK % Returns the implicit stiffness field as a
            %   surfaceField. The function computePiolaField should have
            %   been called for an updated impK.
            
            surfaceImpK = obj.materialLaw().surfaceImpK();
        end
        
        function volImpK = volImpK(obj)
            % volImpK % Returns the implicit stiffness field as a volField.
            %   The function computePiolaField should have been called for
            %   an updated impK.
            
            volImpK = obj.materialLaw().volImpK();
        end
    end
end

