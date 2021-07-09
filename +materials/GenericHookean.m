classdef GenericHookean < materials.Hookean
    %GenericHookean Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Switch for plane stress condition
        planeStress_
        
        % Density
        density_
        
        % Lame's first parameter
        lambda_;
        
        % Shear modulus (also known as G)
        mu_
    end
    
    methods
        function obj = GenericHookean(planeStress, params)
            %GenericHookean Construct an instance of this class
            %   Detailed explanation goes here
            
            [ lambda, mu ] = ...
                computeLameConstants...
                (...
                    struct...
                    (...
                        'E', params.E,...
                        'Nu', params.Nu,...
                        'planeStress', planeStress ...
                    )...
                );
            
            obj.planeStress_ = planeStress;
            obj.density_ = params.density;
            obj.lambda_ = lambda;
            obj.mu_ = mu;
        end
        
        %% BEG - Hookean interface implementation
        function planeStress = planeStress(obj)
            % planeStress Returns the switch for plane stress condition.
            
            planeStress = obj.planeStress_;
        end
        
        function obj = setPlaneStress(obj, planeStress)
            % setPlaneStress Sets the switch for plane stress condition.
            
             obj.planeStress_ = planeStress;
        end
        
        function density = density(obj)
            % density Density
           
            density = obj.density_;
        end
        
        function obj = setDensity(obj, density)
            % setDensity Set the Density
            
            obj.density_ = density;
        end
        
        function lambda = lambda(obj)
            % lambda The Lame's first parameter
            
            lambda = obj.lambda_;
        end
        
        function obj = setLambda(obj, lambda)
            %lambda Set the Lame's first parameter
            
            obj.lambda_ = lambda;
        end
        
        function mu = mu(obj)
            % mu is the Shear Modulus
            
            mu = obj.mu_;
        end
        
        function obj = setMu(obj, mu)
            % setMu Set the Shear Modulus
            
            obj.mu_ = mu;
        end
        %% END - Hookean interface implementation
    end
end

