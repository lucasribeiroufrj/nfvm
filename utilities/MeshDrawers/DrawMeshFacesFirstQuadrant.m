classdef DrawMeshFacesFirstQuadrant < MeshDrawer & handle
    %meshDrawer Summary of this class goes here
    %   Detailed explanation goes here
    %
    % TODO:
    % - Only works with mesh having faces with the same number of vertices.
    %
    
    properties
    end
    
    methods
        function obj = DrawMeshFacesFirstQuadrant(mesh, figHandle)
            %DrawMeshFacesAtFirstQuadrant Construct an instance of this 
            % class detailed explanation goes here
            
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
            faces = zeros(mesh.numberOfFaces, nVert);
            vertices = zeros(mesh.numberOfNodes, 3);
            %%
            
            %% BEG - Collect faces which satisfies a criteria.
            %
            faceCounter = 0;
            vertices(1:end) = NaN; %Flag vert. not to me included.

            for iFace = 1:mesh.numberOfFaces

                nodes = mesh.faces(iFace).iNodes;
                satisfiedCriteria = false;

                for j = 1:nVert

                    iVert = nodes(j);
                    c = mesh.nodes(iVert).centroid;

                    if c(1) < 0 || c(2) < 0 || c(3) < 0
                        continue
                    else
                        satisfiedCriteria = true;
                        vertices(iVert,:) = c;
                    end
                end

                if satisfiedCriteria
                    faceCounter = faceCounter + 1;
                    faces(faceCounter,:) = mesh.faces(iFace).iNodes;
                end
            end

            if faceCounter == 0
                disp('No faces in the first quadrant.')
            end

            faces = faces(1:faceCounter,:);
            %% END
        end
    end
end

