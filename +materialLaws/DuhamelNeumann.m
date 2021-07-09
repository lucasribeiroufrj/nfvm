classdef DuhamelNeumann < materialLaws.HookeanElastic
    %DuhamelNeumann Law described in [eq. 7, Demirdzic - 1994].
    %   Detailed explanation goes here
    %
    %   This is a linear model, so Sigma (second piola) and Piola (first)
    % are equal.
    %
    % Ref.:
    % [Demirdzic - 1994] - Numerical method for coupled fluid flow, Heat 
    % transfer and stress analysis using unstructured moving meshes with 
    % cells of arbitrary topology
    
    properties
        alpha_;
        deltaTemp_;
    end
    
    methods
        function obj = ...
                DuhamelNeumann...
                (...
                    mesh,...
                    material,...
                    alpha,...
                    deltaTemp...
                )
            
            obj@materialLaws.HookeanElastic...
            (...
                mesh,...
                material...
            );
        
            obj.alpha_ = alpha;
            obj.deltaTemp_ = deltaTemp;
        end
        
        %% BEG - materialLaws.SegregatedBaseImpl methods implementation
        %
        function Sigma = computeSigmaTensorField(obj, gradU)
            % Computes the second piola stress given the right-Cauchy-Green
            %   strain tansor (as a tensorField). The equation is described 
            %   in [6.28, Bonet - 2008].
            %
            %   computeSigmaTensorField(gradU) returns a tensorField. The
            %       argument gradU must be tensorField too.
            
            mu = obj.material().mu();
            lambda = obj.material().lambda();
            alpha = obj.alpha_;
            deltaTemp = obj.deltaTemp_;
            
            HookeanStress = ...
                obj.computeSigmaTensorField@materialLaws.HookeanElastic...
                (...
                    gradU...
                );
            
            I = gradU.identity();
            Sigma = HookeanStress - (2*mu + 3*lambda)*alpha*deltaTemp*I;
        end
        %% END - materialLaws.SegregatedBaseImpl methods implementation
    end
end

