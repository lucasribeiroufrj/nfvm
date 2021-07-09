classdef HookeanElastic ...
        < materialLaws.BlockCoupledMaterialLaw...
        & materialLaws.NonLinearBlockCoupledBaseImpl ...
        & materialLaws.SegregatedBaseImpl ...
        & handle
    %HookMaterial Hook's law [eq. 7, Demirdzic - 1994].
    %   This material law can be used with BC, Segregated and NLBC.
    %
    % Ref.:
    % [Demirdzic - 1994] - Numerical method for coupled fluid flow, Heat 
    % transfer and stress analysis using unstructured moving meshes with 
    % cells of arbitrary topology - Demirdzic - 1994
    
    properties
        
    end
    
    methods
        
        function obj = HookeanElastic(mesh, material)
            % HookeanElastic creates a HookeanElastic material law.
            %   mesh : must be fvMesh-compatible.
            %   material : must be an HookeanElastic-compatible material.
        
            obj@materialLaws.SegregatedBaseImpl(mesh, material); 
        end
        
        %% BEG - BlockCoupledMaterialLaw methods implementation
        %
        function mu = mu(obj)
            % mu Shear modulus (also known as G).
            
            mu = obj.material().mu();
        end
        
        function lambda = lambda(obj)
            % lambda Lame's first parameter.
           
            lambda = obj.material().lambda();
        end
        %% END - BlockCoupledMaterialLaw methods implementation
        
        
        %% BEG - materialLaws.SegregatedBaseImpl methods implementation
        %        
        function Sigma = computeSigmaTensorField(obj, gradU)
            % Computes the second piola stress for a given displacement
            %   gradient field (as tensorField). The equation is described 
            %   in [6.28, Bonet - 2008].
            %
            %   computeSigmaTensorField(gradU) returns a tensorField. The
            %       argument gradU must be tensorField too.
            
            I = gradU.identity();
            
            mu = obj.material().mu();
            lambda = obj.material().lambda();

            e = (0.5) * (gradU + gradU.');
            
            tr_e = e.trace();
            Sigma = 2*mu*e + lambda*tr_e*I;
        end
        %% END - materialLaws.SegregatedBaseImpl methods implementation
        
        
        %% BEG - NonLinearBlockCoupledBaseImpl interface implementation
        %
        function Td = computeTdTensorField(obj, F, n, t)
            % computeTdTensorField returns Td = FC : nt
            %                                    (2,3)
            % where C is elasticity tensor, n is the normal face and t is
            % any vector. The inputs must be tensor, vector and vector
            % fields respectively.
            
            mu = obj.material().mu();
            lambda = obj.material().lambda();
            Td = lambda*((F*n)&t) + mu*((t^n)*F + ((F*t)&n));
        end
        %% END - NonLinearBlockCoupledBaseImpl methods implementation
    end
end

