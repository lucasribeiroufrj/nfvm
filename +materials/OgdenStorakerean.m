classdef (Abstract) OgdenStorakerean < materials.Material
    %Hookean Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods (Abstract)

        alphas = alphas(obj)
        betas = betas(obj)
        mus = mus(obj)
    end
end

