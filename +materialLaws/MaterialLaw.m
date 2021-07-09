classdef (Abstract) MaterialLaw < handle
    %MaterialLaw Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Abstract)
        
        % Computes the first-piola stress tensor for a given displacement 
        % gradient field. The gradU parameter can be a volume or a surface
        % fields.
        Piola = computePiolaField(obj, gradU)
        
        % Returns the density
        density = density(obj)
        
        % A boolean to indicate whether the stress is planar or not.
        planeStress = planeStress(obj)
    end
    
end

