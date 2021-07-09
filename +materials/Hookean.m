classdef Hookean < materials.Material
    %Hookean Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods (Abstract)
        
        lambda = lambda(obj)        
        obj = setLambda(obj, lambda)        
        mu = mu(obj)        
        obj = setMu(obj, mu)
    end
end

