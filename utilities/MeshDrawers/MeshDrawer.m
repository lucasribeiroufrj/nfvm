classdef MeshDrawer < handle
    %meshDrawer Summary of this class goes here
    %   Detailed explanation goes here
    %
    % TODO:
    % - Only works with mesh having faces with the same number of vertices.
    %
    
    properties
        figHandle_ = []
        mesh_
        faces_
        vertices_
        patch_ = []
        vertData_;
    end
    
    methods
        function obj = MeshDrawer(mesh, figHandle)
            %meshDrawer Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.mesh_ = mesh;
            
            if nargin == 2
                obj.figHandle_ = figHandle;
            end
            
            [obj.faces_, obj.vertices_] = obj.collectFaces();
            
            obj.initialization();
        end
        
        function vertices = displacedVertices(obj, volPhi, volUGrad)
            
            vertices = obj.vertices_;
            
            phi = volPhi.internalField().data;
            gradPhi = volUGrad.internalField().data;
            spaceMesh = obj.mesh_.spaceMesh();
            
            for iVert = 1:size(obj.vertices_,1) 

                % Satisfied criteria?
                if ~isnan(obj.vertices_(iVert,1))

                   iElems = spaceMesh.nodes(iVert).iElements(:);
                   totalU = [0; 0; 0];
                   
                   for idx = 1:length(iElems)
                       
                       iElem = iElems(idx);
                       
                       invDist = obj.vertData_(iVert).invDist(iElem);
                       d_vc = obj.vertData_(iVert).d_vc(:,iElem);
                       
                       interpU = phi(:, iElem) + gradPhi(:,:,iElem)*d_vc;
                       totalU = totalU + invDist*interpU;
                   end
                   
                   totalU = obj.vertData_(iVert).weight*totalU;
                   vertices(iVert,:) = vertices(iVert,:) + totalU.';
                end
            end
        end
        
        function obj = initialization(obj)
            
            mesh = obj.mesh_;
            vertices = obj.vertices_;
            centroids = mesh.r_C.internalField().data;
            
            nVertices = size(obj.vertices_,1) ;
            
            obj.vertData_(nVertices).d_vc = [];
            
            
            for iVert = 1:nVertices

                % Satisfied criteria?
                if ~isnan(vertices(iVert,1))

                   iElems = mesh.nodes(iVert).iElements(:);
                   factor = 0;
                   
                   obj.vertData_(iVert).d_vc = ...
                       zeros(3,length(iElems));
                   
                   verticePosition = vertices(iVert,:).';
                   
                   for idx = 1:length(iElems)
                       
                       iElem = iElems(idx);
                       d_vc = verticePosition - centroids(:,iElem);
                       invDist = 1 / norm(d_vc);
                       obj.vertData_(iVert).invDist(iElem) = invDist;
                       obj.vertData_(iVert).d_vc(:,iElem) = d_vc;
                       %interpU = phi(:, iElem) + gradPhi(:,:,iElem)*d_vc;
                       %totalU = totalU + invDist*interpU;
                       factor = factor + invDist;
                   end

                   obj.vertData_(iVert).weight = 1/factor;
                end
            end
        end
        
        function run(obj, volPhi, volUGrad)
            % run Draw the mesh (patch)
            %   If only volPhi is passed, then the mesh is drawn without 
            %   changes. But if the volume displacement gradiend volUGrad
            %   is also passed then the mesh will use this field to draw
            % the deformed mesh.
            
            if nargin == 3
                vertices = obj.displacedVertices(volPhi, volUGrad);
            else
                vertices = obj.vertices_;
            end
            
            if isempty(obj.figHandle_)
                obj.figHandle_ = figure(999);
                clf(obj.figHandle_);
                set...
                (...
                    obj.figHandle_,...
                    'NumberTitle', 'off', ...
                    'Name', ...
                    'Mesh'...
                );
            end
            
            % Erase last draw
            if ~isempty(obj.patch_)
                delete(obj.patch_)
            end
            
            figure(obj.figHandle_);
            
            obj.patch_ = ...
                patch...
                (...
                    'Faces',obj.faces_,...
                    'Vertices',vertices,...
                    'FaceColor','green'...
                );
            
            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1]...%,...
                ...%'XLim',[-0.1 1.1],...
                ...%'YLim',[-0.1 1.1],...
                ...%'ZLim',[-0.1 1.1]...
             )
         
            xlabel('x') 
            ylabel('y')
            zlabel('z')

            grid on
            drawnow
        end
    end
    
    methods (Abstract)
        
        [faces, vertices] = collectFaces(obj);
    end
end

