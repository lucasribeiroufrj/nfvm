classdef BoundaryMeshDrawer < MeshDrawer & handle
    %meshDrawer Summary of this class goes here
    %   Detailed explanation goes here
    %
    % TODO:
    % - Only works with mesh having faces with the same number of vertices.
    %
    
    properties
    end
    
    methods
        function obj = BoundaryMeshDrawer(mesh, figHandle)
            %meshDrawer Construct an instance of this class
            %   Detailed explanation goes here
            
            args{1} = mesh;
            
            if nargin == 2
                args{2} = figHandle;
            end
            
            obj@MeshDrawer(args{:});
        end
        
        function [faces, vertices] = collectFaces(obj)
            %collectFaces Collect faces which satisfies a criteria.
            %   Detailed explanation goes here
            
            % Take a reference to the mesh
            mesh = obj.mesh_;
            
            %% Allocate space to store the faces which will be drawn.
            %
            nVert = length(mesh.faces(1).iNodes);
            faces = zeros(mesh.numberOfBFaces, nVert);
            vertices = zeros(mesh.numberOfNodes, 3);
            %%
            
            %% BEG - Collect boundary faces
            %
            vertices(1:end) = NaN; % Flags vertices as internal.

            for i = 1:mesh.numberOfBFaces
                
                faces(i,:) = ...
                    mesh.faces(mesh.numberOfInteriorFaces+i).iNodes;

                for j = 1:nVert

                    idxVert = faces(i,j);
                    vertices(idxVert,:) = mesh.nodes(idxVert).centroid;
                end
            end
            %% END
        end
        
        function obj = initialization2(obj)
            
            mesh = obj.mesh_;
            vertices = obj.vertices_;
            
            nVertices = size(obj.vertices_,1) ;
            
            obj.vertData_(nVertices).d_vc = [];
            
            for iVert = 1:nVertices

                % Satisfied criteria?
                if ~isnan(vertices(iVert,1))

                   iFaces = mesh.nodes(iVert).iFaces(:);
                   factor = 0;
                   
                   obj.vertData_(iVert).d_vc = ...
                       zeros(3,length(iFaces));
                   
                   verticePosition = vertices(iVert,:).';
                   
                   for idx = 1:length(iFaces)
                       
                       iFace = iFaces(idx);
                       d_vc = verticePosition - mesh.faces(iFace).centroid;
                       invDist = 1 / norm(d_vc)^2;
                       obj.vertData_(iVert).invDist(iFace) = invDist;
                       obj.vertData_(iVert).d_vc(:,iFace) = d_vc;
                       factor = factor + invDist;
                   end

                   obj.vertData_(iVert).weight = 1/factor;
                end
            end
        end
        
        function vertices = displacedVertices2(obj, volPhi, ~)
            
            vertices = obj.vertices_;
            mesh = obj.mesh_.spaceMesh();
            
            for iVert = 1:size(obj.vertices_,1) 

                % Satisfied criteria?
                if ~isnan(obj.vertices_(iVert,1))
                    
                    iFaces = mesh.nodes(iVert).iFaces(:);
                    totalU = [0; 0; 0];
                    
                    for idx = 1:length(iFaces)
                       
                       iFace = iFaces(idx);
                       face = mesh.faces(iFace);
                       iBoundary = face.patchIndex;
                       
                       if iBoundary == 0; continue; end
                       
                       boundary = volPhi.boundaryField(iBoundary);
                       u = boundary.field(iFace - boundary.startFace + 1);
                       invDist = obj.vertData_(iVert).invDist(iFace);
                       totalU = totalU + invDist*u;
                    end
                   
                   totalU = obj.vertData_(iVert).weight*totalU;
                   vertices(iVert,:) = vertices(iVert,:) + totalU.';
                end
            end
        end
    end
end

