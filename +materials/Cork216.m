classdef Cork216 < materials.GenericHookean & materials.OgdenStorakerean
    %Cork216 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
                
        % OgdenStorakers paramters
        alphas_
        betas_
        mus_
    end
    
    methods
        function obj = Cork216(planeStress)
            %Cork216 Construct an instance of this class
            %   Detailed explanation goes here
            
            params.E = 20e+6;
            params.Nu = 0.3;
            params.density = 216;
            
            obj = obj@materials.GenericHookean(planeStress, params);
            
            % OgdenStorakerean parameters
            obj.alphas_ = [19.35 28.58 19.9];
            obj.betas_ = [0.5423 2.0905 0.3967];
            obj.mus_ = [1.982e+6 0.00015478e+6 1.1045e+6];
        end
        
        
        %% BEG - OgdenStorakerean interface implementation
        function alphas = alphas(obj)
            
            alphas = obj.alphas_;
        end
        
        function betas = betas(obj)
            
            betas = obj.betas_;
        end
        
        function mus = mus(obj)
            
            mus = obj.mus_;
        end
        %% END - OgdenStorakerean interface implementation
    end
end

