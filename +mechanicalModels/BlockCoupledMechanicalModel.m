classdef BlockCoupledMechanicalModel...
        < mechanicalModels.MechanicalModel & handle
    %BlockCoupledMechanicalModel Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)

        % The MaterialLaw. Might be Hookean, OgdenStorakerean, etc.
        materialLaw_;
    end
    
    methods
        function obj = BlockCoupledMechanicalModel(materialLaw)
            % BlockCoupledMechanicalModel Construct an mechanical model.
            %   This mechanical model is compatible with the BlockCoupled
            %   solid solver.
            %   materialLaw : might be Hookean or any other 
            %   BlockCoupledMaterialLaw-derived material law.
            
            obj.materialLaw_ = materialLaw;
        end
        
        % Beg - mechanicalModels.MechanicalModel interface implementation
        %
        function Piola = computePiolaField(obj, gradUField)
            % computePiolaField return the first-piola stress tensor a 
            %   given displacement gradient field. The gradUField parameter
            %   can be a volume or a surface fields.
            
            Piola = obj.materialLaw().computePiolaField(gradUField);
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
        
        function mu = mu(obj)
            % mu Shear modulus (also known as G).
            
            mu = obj.materialLaw().mu();
        end
        
        function lambda = lambda(obj)
            % lambda Lame's first parameter.
           
            lambda = obj.materialLaw().lambda();
        end
    end
end

