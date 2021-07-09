classdef NeoHookean ...
        < materialLaws.NonLinearBlockCoupledBaseImpl ...
        & materialLaws.SegregatedBaseImpl ...
        & handle
    %HookMaterial Neo-Hookean's law [Bonet - 2008].
    %   This material law can be used with Segregated and NLBC.
    %
    % Ref.:
    % [Bonet - 2008] - Nonlinear continuum mechanics for finite element
    % analysis.
    
    properties
       
    end
    
    methods
        
        function obj = NeoHookean(mesh, material)
            % NeoHookean creates a NeoHookean material law.
            %   mesh : must be fvMesh-compatible.
            %   material : must be an NeoHookean-compatible material.
            
            obj@materialLaws.SegregatedBaseImpl(mesh, material);
        end
        
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
            F = I + gradU;
            C = F.'*F;
            
            lnJ = C.det().sqrt_().log_();
            invC = C.inv();
            mu = obj.material().mu();
            lambda = obj.material().lambda();
            Sigma = mu*(I-invC) + (lambda*lnJ)*invC;
        end
        %% END - materialLaws.SegregatedBaseImpl methods implementation
        
        
        %% BEG - NonLinearBlockCoupledBaseImpl methods implementation
        %
        function Td = computeTdTensorField(obj, F, n, t)
            % computeTdTensorField returns Td = FC : nt
            %                                    (2,3)
            % where C is elasticity tensor, n is the normal face and t is
            % any vector. The inputs must be tensor, vector and vector
            % fields respectively.
            
            C = F.'*F;
            lnJ = fields.log(fields.sqrt(fields.det(C)));
            invC = C.inv();
            
            mu = obj.material().mu();
            lambda = obj.material().lambda();
            
            one = lnJ.ones();
            A = F*invC;
            b = invC*n;
            Td = ((A*(lambda*n)) & (invC*t)) + ...
                (mu*one - lambda*lnJ)*(((A*t) & b) + (b^t)*A);
        end
        %% END - NonLinearBlockCoupledBaseImpl methods implementation
    end
end
