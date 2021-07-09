classdef (Abstract) SolidModel < handle
    %SolidModel Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Finite volume mesh:
        fvMesh_;
    end
    
    methods
        function obj = SolidModel(fvMesh)
            %SolidModel Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.fvMesh_ = fvMesh;
        end
        
        %% BEG - Data member access
        function fvMesh = fvMesh(obj)
           
            fvMesh = obj.fvMesh_;
        end
        
        function spaceMesh = mesh(obj)
           
            spaceMesh = obj.fvMesh_.spaceMesh();
        end
        %% END
    end
    
    methods (Abstract)
       
        % Evolve the solid model.
        obj = evolve(obj)
        
        % Returns the first piola stress volume tensor.
        volPiola = volPiola(obj)
        
        % Returns the first piola stress surface tensor.
        surfacePiola = surfacePiola(obj)
    end
end

