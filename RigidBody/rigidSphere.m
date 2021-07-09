classdef rigidSphere
    %rigiSphere Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Old position
        position_o_
        
        % Current position
        position_
        
        % Velocity function of time.
        velocity_
        
        % Radius of the sphere
        radii_
    end
    
    methods
        function obj = ...
                rigiSphere...
                (...
                    position_o,...
                    position,...
                    velocity,...
                    radii...
                )
            %rigiSphere Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.position_o_= position_o;
            obj.position_ = position;
            obj.velocity_ = velocity;
            obj.radii_ = radii;
        end
        
        function obj = move(obj, deltaT)
            %move Integrate the velocity in time.
            %   New position = old position + velocity*deltaT
            
            obj.position_ = obj.position_o_ + obj.velocity_*deltaT;
        end
    end
end

