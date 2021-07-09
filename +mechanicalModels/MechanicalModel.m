classdef (Abstract) MechanicalModel < handle
    %MechanicalModel Summary of this class goes here
    %   Detailed explanation goes here

    methods (Abstract)
        
        % Computes the first-piola stress tensor for a given displacement 
        % gradient field. The gradU parameter can be a volume or a surface
        % fields.
        Piola = computePiolaField(obj, gradU)
        
        % Returns the material law
        materialLaw = materialLaw(obj)
        
        % Returns the density
        density = density(obj)
    end
end

