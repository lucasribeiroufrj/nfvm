classdef rigidSphereContact
    %rigidSphereContact Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Radius of the sphere
        radii_;
        
        % Initial position of Sphere's center.
        x0_;
        
        % Impact duration
        tf_;
        
        % Max acceleration
        acc_max_ = -225;
        
        % Polinomial's coefficients
        a_;
        b_;
        
        % Initial velocity;
        v0_;
        
        % Acceleration function
        acceleration_;
        
        % Time interval
        t_;
        
        % Position impact
        position_;
    end
    
    methods
        
        function obj = ...
            rigidSphereContact...
            (...
                radii,...
                initialPosition,...
                initialVelocity,...
                endTime...
            )
            %rigidSphereContact Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.radii_ = radii;
            obj.x0_ = initialPosition; %radii + heightBrick;
            obj.v0_ = initialVelocity;
            
            obj.tf_ = endTime;
            obj.a_ = obj.acc_max_/((obj.tf_/2)^2);
            obj.b_ = -obj.tf_*obj.a_;
            
            obj.t_ = (0:obj.tf_/18:obj.tf_);
            
            a_ = obj.a_;
            v0_ = obj.v0_;
            x0_ = obj.x0_;
            b_ = obj.b_;
            obj.acceleration_ = @(t) a_*t.^2 +b_*t;
            obj.position_ = @(t,v0, x0) a_/12*t.^4 + b_/6*t.^3 + v0_*t + x0_;
        end
        
        function position = position(obj, time)
            
            position = obj.position_(time);
        end
        
        function plotAcceleration(obj)
            
            figure(1);
            plot(obj.t_,obj.acceleration_(obj.t_));
            xlabel('Time [s]');
            ylabel('Acceleration [g]');
        end
        
        function plotPosition(obj)
            
            figure(2)
            clf

            plot...
            (...
                obj.t_,...
                obj.position_...
                (...
                    obj.t_, obj.v0_, obj.x0_...
                )...
            );
            xlabel('Time [s]');
            ylabel('Sphere''s center position [m]');
            hold on
            plot...
            (...
                obj.t_,...
                obj.position_(obj.t_, obj.v0_, obj.x0_)-obj.radii_...
            );

            % Minimum height brick 34.23/2*1e-3
            top = @(t) t*0 + obj.x0_ - obj.radii_;
            plot(obj.t_,top(obj.t_));
            axis([0 9e-3 0 0.09])
        end
        
        function radii = radii(obj)
            
            radii = obj.radii_;
        end
        
        function obj = plotBoth(obj)
            
            obj.plotAcceleration();
            obj.plotPosition();
        end
    end
end

