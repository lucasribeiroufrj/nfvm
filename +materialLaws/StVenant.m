classdef StVenant ...
        < materialLaws.NonLinearBlockCoupledBaseImpl ...
        & materialLaws.SegregatedBaseImpl ...
        & handle
    %StVenant material law [Ciarlet - 1988].
    %   This is the simplest hyperelastic model. It should not be expected
    %   to behave well for large strains, since it is only modeled
    %   according to the behaviour for "small" values of of the 
    %   Green-Lagrange strain tensor ||E|| [p.159, Ciarlet - 1988]. This is
    %   why such law is often referred to as "large displacement-small
    %   strain" model. It be expected to perform better than  the
    %   linearized models that ate so often used [p.132, Ciarlet - 1988].
    %
    %   This material law can be used with Segregated and NLBC.
    %
    % Ref.:
    % [Ciarlet - 1988] - Mathematical Elasticity - Volume I - 
    % Three-Dimensional Elasticity - Ciarlet - 1988
    
    properties

    end
    
    methods
        
        function obj = StVenant(mesh, material)
            % StVenant creates a StVenant material law.
            %   mesh : must be fvMesh-compatible.
            %   material : must be an StVenant-compatible material.
            
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
            
            mu = obj.material().mu();
            lambda = obj.material().lambda();
            
            E = 0.5*(C - I);
            
            trE = E.trace();
            Sigma = 2*mu*E + lambda*trE*I;
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
            
            mu = obj.material().mu();
            lambda = obj.material().lambda();
            Td = lambda*((F*n)&t) + mu*((t^n)*F + ((F*t)&n));
        end
        %% END - NonLinearBlockCoupledBaseImpl methods implementation
    end 
end

