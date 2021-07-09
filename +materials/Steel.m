classdef Steel < materials.GenericHookean
    %Steel Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = Steel(planeStress)
            %Steel Construct an instance of this class
            %   Detailed explanation goes here
            
            params.E = 2e+11;
            params.Nu = 0.3;
            params.density = 1e+4;
            
            obj = obj@materials.GenericHookean(planeStress, params);
        end
    end
end

