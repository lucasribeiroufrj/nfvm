classdef (Abstract) Material
    %Material Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods (Abstract)
        
        planeStress = planeStress(obj)
        obj = setPlaneStress(obj, planeStress)
        density = density(obj)
        obj = setDensity(obj, density)
    end
end

