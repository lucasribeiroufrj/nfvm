classdef Rubber < materials.GenericHookean
    %Rubber Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods
        function obj = Rubber(planeStress)
            %Rubber Construct an instance of this class
            %   Detailed explanation goes here
            
            params.E = 3e+7;
            params.Nu = 0.3;
            params.density = 1e+3;
            
            obj = obj@materials.GenericHookean(planeStress, params);
        end
    end
end

