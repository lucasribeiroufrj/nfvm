classdef LeastSquareVolPointInterpolation < handle
    %LeastSquareVolPointInterpolation Inspiration came from solids4foam.
    %   Inspired by solids4foam's newLeastSquaresVolPointInterpolation.
    %   Every mesh point has:
    %    - a local origin;
    %    - a weighting factor;
    %    - an inverse of the least square linear interpolation matrix;
    %    - a mirror plane transformation.
    
    properties
        
        % Finite volume mesh.
        fvMesh_;
        
        % local origins.
        origins_ = [];
        
        % weighting factors.
        weights_ = [];
        
        % Precomputed coefficients.
        coeffs_ = [];
        
        % inverses of the least square linear interpolation matrices.
        invLsMatrices_ = [];
        
        % Mirror plane normals and transformation tensors
        mirrorPlaneTransformations_ = [];
    end
    
    methods
        
        function obj = LeastSquareVolPointInterpolation(fvMesh)
            
            obj.fvMesh_ = fvMesh;
        end

        function origins = origins(obj)
            % origins() returns an array of local origins.
            
            if isempty(obj.origins_)
                
                obj.origins_ = ...
                    LeastSquareVolPointInterpolation.makeOrigins...
                    (...
                        obj.fvMesh_.spaceMesh(),...
                        obj.mirrorPlaneTransformations(),...
                        obj.weights()...
                    );
            end
            
            origins = obj.origins_;
        end
        
        function weights = weights(obj)
           % weights() returns an array of weighting factors.
            
           if isempty(obj.weights_)
              
               obj.weights_ = ...
                    LeastSquareVolPointInterpolation.makeWeights...
                    (...
                        obj.fvMesh_.spaceMesh(),...
                        obj.mirrorPlaneTransformations()...
                    );
           end
           
           weights = obj.weights_;
        end
        
        function invLsMatrices = invLsMatrices(obj)
           % invLsMatrices() returns an array of inverses of the least
           %    square linear interpolation matrices.
            
           if isempty(obj.invLsMatrices_)
              
               obj.invLsMatrices_ = ...
                    LeastSquareVolPointInterpolation.makeInvLsMatrices...
                    (...
                        obj.fvMesh_.spaceMesh(),...
                        obj.mirrorPlaneTransformations(),...
                        obj.weights(),...
                        obj.origins()...
                    );
           end
           
            invLsMatrices = obj.invLsMatrices_;
        end
        
        function mirrorPlaneTransformations = ...
                mirrorPlaneTransformations(obj)
            % mirrorPlaneTransformations() returns an array of the mirror
            %    plane normals and transformation tensors
            
            if isempty(obj.mirrorPlaneTransformations_)
                
                obj.mirrorPlaneTransformations_ = ...
                    LeastSquareVolPointInterpolation...
                    .makeMirrorPlaneTransformations...
                    (...
                        obj.fvMesh_.spaceMesh()...
                    );
            end
            
            mirrorPlaneTransformations = obj.mirrorPlaneTransformations_;
        end
    end
    
    methods (Static)
        
        function origins = ...
                makeOrigins(mesh, mirrorPlaceTransformations, weights)
            % makeOrigins() returns an array of local origins.
            
            nNodes = mesh.numberOfNodes();
            origins = zeros(3,nNodes);
            
            for n = 1:nNodes
               
                node = mesh.nodes(n);
                iElements = node.iElements;
                iBFaces = node.iBFaces;
                nElements = length(iElements);
                nBFaces = length(iBFaces);
                nAllPoints = nElements + nBFaces;
                
                allPoints = zeros(3, nAllPoints);
                
                for p = 1:nElements
                    
                    centroid = mesh.elements(iElements(p)).centroid;
                    allPoints(:,p) = centroid;
                end
                
                for f = 1:nBFaces
                    
                    centroid = mesh.faces(iBFaces(f)).centroid;
                    allPoints(:,p + f) = centroid;
                end
                
                allMirrorPoints = [];
                nf = mirrorPlaceTransformations(n).normal;
                if norm(nf) > eps
                    
                    allMirrorPoints = zeros(3, nAllPoints);
                    I = eye(3,3);
                    
                    for m = 1:nAllPoints
                       
                        centroid = node.centroid;
                        delta = allPoints(:,m) - centroid;
                        allMirrorPoints(:,m) = ...
                            centroid + (I-(2*nf)*nf.')*(delta);
                    end
                end
                
                normFactor = 0;
                
                for p = 1:nAllPoints
                    
                    W = weights(n).points(p)^2;
                    normFactor = normFactor + W;
                    origins(:,n) = origins(:,n) + W*allPoints(:,p);
                end

                for m = 1:size(allMirrorPoints,2)

                    W = weights(n).points(p + m)^2;
                    normFactor = normFactor + W;
                    origins(:,n) = origins(:,n) + W*allMirrorPoints(:,m);
                end
                
                origins(:,n) = origins(:,n)/normFactor;
            end
        end
        
        function weights = makeWeights(mesh, mirrorPlaneTransformations)
            % makeWeights() returns an array of the mirror plane normals
            %    and transformation tensors
            
            nNodes = mesh.numberOfNodes;
            weights(nNodes) = struct('points', []);
            
            for n = 1:nNodes
                
                node = mesh.nodes(n);
                iElements = node.iElements;
                iBFaces = node.iBFaces;
                nElements = length(iElements);
                nBFaces = length(iBFaces);
                nAllPoints = nElements + nBFaces;
                
                allPoints = zeros(3, nAllPoints);

                for p = 1:nElements

                    centroid = mesh.elements(iElements(p)).centroid;
                    allPoints(:,p) = centroid;
                end

                for f = 1:nBFaces

                    centroid = mesh.faces(iBFaces(f)).centroid;
                    allPoints(:,p + f) = centroid;
                end
                
                allMirrorPoints = [];
                nf = mirrorPlaneTransformations(n).normal;
                if norm(nf) > eps
                    
                    allMirrorPoints = zeros(3, nAllPoints);
                    I = eye(3,3);
                    
                    for m = 1:nAllPoints
                       
                        centroid = node.centroid;
                        delta = allPoints(:,m) - centroid;
                        allMirrorPoints(:,m) = ...
                            centroid + (I-(2*nf)*nf.')*(delta);
                    end
                end
                
                nAllMirrorPoints = size(allMirrorPoints,2);
                total = nAllMirrorPoints + nElements + nBFaces;
                
                % From newLeastSquaresVolPointInterpolation.C@solids4foam:
                % "philipc: force weights to 1.0: required for block 
                % coupled solver also arguably for stable for segregated"
                weights(n).points = ones(1, total);
            end
        end
        
        function invLsMatrices = ...
                makeInvLsMatrices...
                (...
                    mesh,...
                    mirrorPlaceTransformations,...
                    weights,...
                    origins...
                )
           % makeInvLsMatrices() returns an array of inverses of the least
           %    square linear interpolation matrices.
            
            nNodes = mesh.numberOfNodes;
            invLsMatrices(nNodes) = struct('value', []);
            
            for n = 1:nNodes
                
                node = mesh.nodes(n);
                iElements = node.iElements;
                iBFaces = node.iBFaces;
                nElements = length(iElements);
                nBFaces = length(iBFaces);
                nAllPoints = nElements + nBFaces;
                
                allPoints = zeros(3, nAllPoints);
                
                for p = 1:nElements
                    
                    centroid = mesh.elements(iElements(p)).centroid;
                    allPoints(:,p) = centroid;
                end
                
                for f = 1:nBFaces
                    
                    centroid = mesh.faces(iBFaces(f)).centroid;
                    allPoints(:,p + f) = centroid;
                end
                
                allMirrorPoints = [];
                nf = mirrorPlaceTransformations(n).normal;
                if norm(nf) > eps
                    
                    allMirrorPoints = zeros(3, nAllPoints);
                    I = eye(3,3);
                    
                    for m = 1:nAllPoints
                       
                        centroid = node.centroid;
                        delta = allPoints(:,m) - centroid;
                        allMirrorPoints(:,m) = ...
                            centroid + (I-(2*nf)*nf.')*(delta);
                    end
                end
                
                nAllMirrorPoints = size(allMirrorPoints,2);
                M = zeros(nAllPoints + nAllMirrorPoints,3);
                invLsMatrices(n).value = M.';
                
                for p = 1:nAllPoints
                    
                    V = allPoints(:,p) - origins(:,n);
                    M(p,:) = V.';
                end
                
                for mP = 1:nAllMirrorPoints
                    
                    V = allMirrorPoints(:,mP) - origins(:,n);
                    M(p + mP, :) = V.';
                end
                
                for i = 1:size(M,1)
                    
                    M(i,:) = M(i,:) * weights(n).points(i);
                end
                
                %invLsM = inv(M.'*M);   % Option 1
                lsM = M.'*M;            % Option 2
                
                for p = 1:length(weights(n).points)
                    
                    M(p,:) = M(p,:)*weights(n).points(p);
                end
                
                %invLsMatrices(n).value = invLsM*M.';   % Option 1
                invLsMatrices(n).value = lsM\M.';       % Option 2
            end
        end
        
        function mPlaneTransf = makeMirrorPlaneTransformations(mesh)
            % makeMirrorPlaneTransformations() returns an array of the
            %    mirror plane normals and transformation tensors
            
            nNodes = mesh.numberOfNodes;
            mPlaneTransf(1:nNodes) = ...
                struct('normal', zeros(3,1), 'tensor', zeros(3,3));
            
            for b = 1:mesh.numberOfBoundaries;
                
                boundary = mesh.boundaries(b);
                
                if strcmp(boundary.type, 'empty')
                    
                    pointNormals = boundary.pointNormals;
                    
                    for n = 1:length(boundary.iNodes)
                        
                        iNode = boundary.iNodes(n);
                        
                        mPlaneTransf(iNode).normal = pointNormals(:,n);
                        mPlaneTransf(iNode).tensor = eye(3,3);
                    end
                end
            end
        end
    end
end
