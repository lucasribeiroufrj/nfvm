classdef BlockCoupled < SolidModelSolver & handle
    %BlockCoupled Summary of this class goes here
    %   This class implements the algorithm described in [Cardiff - 2016]
    %
    % Ref.:
    % [Cardiff - 2016] - A block-coupled Finite Volume methodology for
    % linear elasticity and unstructured meshes - Cardiff.
    % [Clifford - 2011] - Block-Coupled Simulations Using OpenFOAM - 
    % Clifford
    
    properties
        
        % To compute tangential discretization.
        lsVolPointInterp_;
        
        % The number of degrees of freedom of the system.
        degreesOfFreedom_;
    end
    
    methods
        function obj = ...
                BlockCoupled(mechanicalModel, problem)
            %BlockCoupled Construct an instance of this class
            %   Detailed explanation goes here
            
            if ~isa...
                    (...
                        mechanicalModel.materialLaw(),...
                        'materialLaws.HookeanElastic'...
                    )
                error('Only HookeanElastic is supported');
            end
            
            obj@SolidModelSolver(mechanicalModel, problem)

            obj.lsVolPointInterp_ = ...
                LeastSquareVolPointInterpolation(obj.fvMesh());
            obj.setDegreesOfFreedom();
        end
        
        %% BEG - Data member access     
        function obj = setDegreesOfFreedom(obj)
            
            mesh = obj.fvMesh().spaceMesh();
            obj.degreesOfFreedom_ = 3*(...
                      mesh.numberOfElements ...
                    + mesh.numberOfBElements ...
                    - numberOfFacesInEmptyBoundaries(mesh)...
                    );
        end
        
        function volPiola = volPiola(obj)
            %volPiola First piola stress tensor.
            %   volPiola() returns the first piola stress volume tensor.
            
            error('To be implemented');
        end
        
        function surfacePiola = surfacePiola(obj)
            %surfacePiola First piola stress tensor.
            %   surfacePiola() returns the first piola stress surface
            %   tensor.
            
            error('To be implemented');
        end
        %% END
        
        % Solve the momentum equation.
        function evolve(obj)
            
            mesh = obj.fvMesh().spaceMesh();
            mu = obj.mechanicalModel.mu();
            lambda = obj.mechanicalModel.lambda();
            
            % Store previews iteration. Needed for residual computation.
            volU_prevIter = obj.volU();
            
            for iCorrector=1:obj.config().nCorrectors
                
                fprintf...
                (...
                    'Iteration %d/%d\n',...
                    iCorrector,...
                    obj.config().nCorrectors...
                );
                
                %% BEG - Compute equation [(6), Cardiff - 2016].
                %
                A = zeros(obj.degreesOfFreedom_,obj.degreesOfFreedom_);
                B = zeros(obj.degreesOfFreedom_,1);
                tic
                % Beg - Algorithm [(1), Cardiff - 2016].
                %
                % Contribution from internal faces.
                for iFace = 1:mesh.numberOfInteriorFaces
                    
                    face = mesh.faces(iFace);
                    
                    coeff = BlockCoupled. ...
                        normalDerivativeCoeff(face, mu, lambda);
                    
                    % beg - Add contribution to matrix diagonal. 
                    i = getTensorialIndex(face.iOwner);
                    A(i, i) = A(i, i) - coeff;
                    
                    j = getTensorialIndex(face.iNeighbour);
                    A(j, j) = A(j, j) - coeff;
                    % end
                    
                    % Add contribution to matrix off-diagonal.
                    A(i, j) = A(i, j) + coeff;
                    A(j, i) = A(j, i) + coeff;
                    
                end % Contribution from internal faces.
                
                % Contribution from boundary faces.
                for iBoundary=1:mesh.numberOfBoundaries
                    
                    boundary = mesh.boundaries(iBoundary);
                    
                    if strcmp(boundary.type,'empty')
                        continue;
                    end

                    for iFace = boundary.startFace:boundary.endFace

                        face = mesh.faces(iFace);
                        
                        coeff = BlockCoupled. ...
                            normalDerivativeCoeff(face, mu, lambda);

                        % Add contribution to matrix diagonal.
                        i = getTensorialIndex(face.iOwner);
                        A(i, i) = A(i, i) - coeff;

                        % Add contribution to matrix off-diagonal.
                        j = getTensorialIndex(face.iBIndex);
                        A(i, j) = A(i, j) + coeff;
                    end
                end  % Contribution from boundary faces.
                %
                % End - Algorithm [(1), Cardiff - 2016].
                
                % Beg - Algorithm [(2), Cardiff - 2016].
                %
                invLsMatrices = obj.lsVolPointInterp_.invLsMatrices();
                weights = obj.lsVolPointInterp_.weights();
                origins = obj.lsVolPointInterp_.origins();
                mPlaneTransf = ...
                    obj.lsVolPointInterp_.mirrorPlaneTransformations();
                
                % For each cell, it computes the tangential derivative.
                for iElement = 1:mesh.numberOfElements
                    
                    iFaces = mesh.elements(iElement).iFaces;
                    l = getTensorialIndex(iElement);
                    
                    for i = 1:length(iFaces)
                        
                        iFace = iFaces(i);
                        face = mesh.faces(iFace);
                        
                        isBoundaryFace = face.patchIndex > 0;

                        if isBoundaryFace
                            
                            boundary = mesh.boundaries(face.patchIndex);
                            
                            if strcmp(boundary.type,'empty')
                                continue;
                            end
                        end
                        
                        n = face.nf;
                            
                        if iElement ~= face.iOwner
                            n = -face.nf;
                        end

                        for edge=face.edges

                            mn = edge.biNormal * n.';
                            mn_t = mn.';
                            L = edge.norm;

                            for j = 1:length(edge.iNodes)
                                
                                %ID for neighbour interpolation point
                                neiID = 1;

                                iNode = edge.iNodes(j);
                                node = mesh.nodes(iNode);
                                dr = node.centroid - origins(:,iNode);
                                iPointElements = ...
                                    [ ...
                                        node.iElements ...
                                        node.iBElements ...
                                    ];
                                invLsMatrix = invLsMatrices(iNode).value;
                                sumW = sum(weights(iNode).points);

                                for k = 1:length(iPointElements)
                                    
                                    M = invLsMatrix(:,neiID);
                                    w = weights(iNode).points(neiID)/sumW;
                                    
                                    coeff = ...
                                        BlockCoupled. ...
                                        tangentialDerivativeCoeff...
                                        (...
                                            L, w, M, dr, mu, lambda, mn,...
                                            mn_t...
                                        );

                                    normNf = ...
                                        norm(mPlaneTransf(iNode).normal);
                                    
                                    if normNf > eps

                                        T = mPlaneTransf(iNode).tensor;
                                        coeff = coeff + T*coeff;
                                    end

                                    iPElement = iPointElements(k);

                                    if iPElement == iElement
                                        % Add contribution to matrix
                                        % diagonal 
                                        A(l, l) = A(l, l) + coeff;
                                    else
                                        % Add contribution to matrix
                                        % off-diagonal 
                                        c = getTensorialIndex(iPElement);
                                        A(l, c) = A(l, c) + coeff;
                                    end

                                    neiID = neiID + 1;
                                end % For each point-cell neighbour
                            end % For each node
                        end % For each edge
                    end % For each cell face
                end % For each cell, it computes the tangential derivative.
                %
                % End - Algorithm [(2), Cardiff - 2016].

                
                % Beg - Contribution from boundary condition
                %
                for iBoundary=1:mesh.numberOfBoundaries
                    
                    boundary = mesh.boundaries(iBoundary);
                    boundaryU = obj.volU().boundaryField(iBoundary);
                    
                    % TODO: <<==== HARD - CODED!!!!
                    % I still don't know to to handle it elegantly.
                    if BlockCoupled. ...
                            isTractionBoundaryCondition(boundaryU)
                    
                        boundaryU = boundaryU.update();
                        
                        externalTraction = ...
                            boundaryU.externalTraction().field();

                        t = 1;
                        for face = mesh.faces(boundary.iFaces)
                            
                            % beg - Normal derivative computation
                            coeff = BlockCoupled. ...
                                normalDerivativeCoeff(face, mu, lambda);
                            
                            % Add contribution to matrix diagonal 
                            l = getTensorialIndex(face.iBIndex);
                            A(l, l) = A(l, l) - coeff;
                            
                            % Add contribution to matrix off-diagonal 
                            c = getTensorialIndex(face.iOwner);
                            A(l, c) = A(l, c) + coeff;
                            % end - Normal derivative computation
                            
                            traction = externalTraction.data(:,t);
                            t = t + 1;
                            sf = norm(face.Sf);
                            B(l) = B(l) - traction*sf;
                            
                            % BEG - Tangential derivative computation
                            %
                            n = face.nf;
                            
                            for edge=face.edges
                                
                                mn = edge.biNormal * n.';
                                mn_t = mn.';
                                L = edge.norm;

                                for node = mesh.nodes(edge.iNodes)
                                    
                                    %ID for neighbour interpolation point
                                    neiID = 1;
                                    iNode = node.index;
                                    dr = node.centroid - origins(:,iNode);
                                    iPointElements = ...
                                        [ ...
                                            node.iElements ...
                                            node.iBElements ...
                                        ];
                                    invLsMatrix = ...
                                        invLsMatrices(iNode).value;
                                    sumW = sum(weights(iNode).points);
                                    
                                    for k = 1:length(iPointElements)

                                        M = invLsMatrix(:,neiID);
                                        wP = weights(iNode).points(neiID);
                                        w = wP / sumW;
                                        
                                        coeff = ...
                                            BlockCoupled. ...
                                            tangentialDerivativeCoeff...
                                            (...
                                                L, w, M, dr, mu, lambda,...
                                                mn, mn_t...
                                            );

                                        Nf = mPlaneTransf(iNode).normal;
                                        if norm(Nf) > eps

                                            T = mPlaneTransf(iNode).tensor;
                                            coeff = coeff + T*coeff;
                                        end

                                        iPElmt = iPointElements(k);

                                        if iPElmt == face.iBIndex
                                            % Add contribution to matrix
                                            % diagonal 
                                            A(l, l) = A(l, l) - coeff;
                                        else
                                            % Add contribution to matrix
                                            % off-diagonal 
                                            c = getTensorialIndex(iPElmt);
                                            A(l, c) = A(l, c) - coeff;
                                        end

                                        neiID = neiID + 1;
                                    end % For each point-cell neighbour
                                end % For each edge node.
                            end
                            % END - Tangential derivative computation
                        end
                        
                    % TODO: <<==== HARD - CODED!!!!
                    % I still don't know to to handle it elegantly.
                    elseif BlockCoupled. ...
                            isDisplacementBoundaryCondition(boundaryU)
                        
                        avgTensor = averageDiagonalTensor(A);
                        diagSign = sum(diag(avgTensor));
                        diagSign = diagSign/norm(diagSign);
                        mag = norm(reshape(avgTensor,1,9));
                        scalarFactor = diagSign*(1/sqrt(3))*mag;
                        T = eye(3,3)*scalarFactor;
                        t = 1;
                        boundaryU = boundaryU.update();
                        
                        for face = mesh.faces(boundary.iFaces)
                            
                            U = boundaryU.field().data(:,t);
                            t = t + 1;
                            l = getTensorialIndex(face.iBIndex);
                            A(l,l) = A(l,l) + T;
                            B(l) = B(l) + scalarFactor*U;
                        end
                        
                    elseif BlockCoupled. ...
                        isSymmetryBoundaryCondition(boundaryU)
                    
                        I = eye(3);
                        
                        for face = mesh.faces(boundary.iFaces)
                            
                            n = face.nf;
                            nn = n * n.';
                            gDiff = norm(face.Sf) / norm(face.CN);
                            T1 = gDiff*(nn - I);
                            coeff_b = -nn + T1;
                            
                            l = getTensorialIndex(face.iBIndex);
                            A(l,l) = A(l,l) + coeff_b;
                            
                            j = getTensorialIndex(face.iOwner);
                            A(l, j) = A(l, j) - T1;
                        end
                        
                    elseif isa(boundaryU, 'emptyPatchVectorField')
                        continue;
                    else
                        ME = MException...
                            (...
                                'BlockCoupled:evolve',...
                                'Invalid boundary patch.'...
                            );
                        throw(ME) 
                    end
                end
                %
                % End - Contribution from boundary condition
                
                %% END - Computes equation [(6), Cardiff - 2016].

                % Add body force
                if obj.problem.hasBodyForce
                    
                    bodyForce = obj.problem.getBodyForce();
                    [r, c] = size(bodyForce.data);
                    dofBodyForce = r*c; 
                    B(1:dofBodyForce) = B(1:dofBodyForce) ...
                        - reshape(bodyForce.data, [dofBodyForce 1]);
                    
                    %obj.plotBodyForce(bodyForce);
                end
                
                if obj.problem.includeAcceleration()
                    
                    [LHS, RHS] = ...
                        BlockCoupled.d2dt2...
                        (...
                            obj.mechanicalModel.density(),...
                            obj.volU_o(),...
                            obj.volU_oo(),...
                            obj.problem()...
                        );

                    for iElem = 1:mesh.numberOfElements
                        
                        l = getTensorialIndex(iElem);
                        dV = obj.mesh.elements(iElem).volume;
                        A(l,l) = A(l,l) - LHS(l,l)*dV;
                        B(l) = B(l) - RHS(l)*dV;
                    end
                end
                
                obj.setVolU(...
                    solveBlockCoupledLinearSystem(obj.volU(), A, B));
                toc
                
                % The 'false' flag is ignore boundary correction.
                obj.correctVolGradU(false); 
                obj.correctSurfaceGradU();

                obj.problem.statistics();
                
                [hasConverged, residual] = ...
                    computeResiduals...
                    (...
                        iCorrector,...
                        obj.volU_,...     % Current iter.
                        volU_prevIter,... % Previews iter.
                        obj.volU_o(),...  % Old timestep.
                        1,... %TODO
                        obj.config()...
                    );

                if hasConverged
                    
                    break;
                end
                
                if isnan(residual)
                    
                    error('NaN detected.')
                end
                
                volU_prevIter = obj.volU_;
            end
            
            % Update fields
            obj.setVolU_oo(obj.volU_o());
            obj.setVolU_o(obj.volU());
        end
        
        function plotBodyForce(obj, bodyForce)
            
            params.mesh = obj.fvMesh;
            maxNorm = max(max(abs(bodyForce.data)));
            bf = obj.fvMesh.r_C()*0;
            bf = bf.setFromRawData...
                (...
                    bodyForce.data/maxNorm*0.1,...
                    1,...
                    bf.nInteriorElements()...
                );
            clf;
            obj.problem.setCamera();
            obj.problem.meshDrawer().run(obj.volU(), obj.volGradU());
            drawnow;
            plotVolVectorField(bf, params);
        end
    end
    
    methods(Static)
        
        function coeff = normalDerivativeCoeff(face, mu, lambda)
            
            delta = norm(face.E);
            d = norm(face.CN);
            n = face.nf;
            nn = n * n.';
            I = eye(3);
            
            coeff = (mu + lambda)*(nn)*(delta / d) + ...
                mu*(I)*(delta / d);
        end
        
        function coeff = ...
                tangentialDerivativeCoeff(L, w, M, dr, mu, lambda, mn,mn_t)
            
            weight = w + (M.' * dr);
            coeff = 0.5*L*weight*(mu*mn + lambda*mn_t);
        end
        
        function [LHS, RHS] = d2dt2(density, U_o, U_oo, problem)
            %fvm_d2dt2 For each finite volume computes the temporal term contribution
            %to the final linear system.
            %
            %   Parameters:
            % - density = the scalar density.
            %
            %   Return:
            % - LHS = The Left-hand-side d2dt2 contribution to the final linear
            % system. 
            % - RHS = The Right-hand-side d2dt2 contribution to the final linear
            % system.
            % 
            %   Obs.: 
            % - The phi_o, phi_oo stands for phi at old, and old old time respectively.

            %% BEG - Initialization
            mesh = U_o.mesh;
            dof = mesh.numberOfElements*3;
            LHS = zeros(dof, dof);

            deltaT = problem.runTime.deltaT();
            deltaT_o = problem.runTime.deltaT_o();
            a = 2*density / ( (deltaT + deltaT_o) * (deltaT * deltaT_o) );
            %% END

            %% BEG - Add contribution from time-dependent term.
            iElements = 1:dof;
            idx = sub2ind(size(LHS), iElements, iElements);

            LHS(idx) = a*deltaT_o;
            
            right = a*( ...
                (deltaT + deltaT_o)*U_o.internalField() ...
                - deltaT*U_oo.internalField() ...
                );
            RHS = reshape(right.data, [dof 1]);
            %% END
        end
    end
end

