classdef NonLinearBlockCoupled < SolidModelSolver & handle
    %NonLinearBlockCoupled Non Linear Block Coupled.
    %   This class implements a non-linear version of the algorithm
    % described in [Cardiff - 2016]. This algorithm implements the Newton-
    % Raphson procedure.
    %
    % Ref.:
    % [Cardiff - 2016] - A block-coupled Finite Volume methodology for
    % linear elasticity and unstructured meshes - Cardiff.
    % [Clifford - 2011] - Block-Coupled Simulations Using OpenFOAM - 
    % Clifford
    
    properties
        
        % Second Piola-Kirchhof
        surfaceSigma_
        
        % To compute tangential discretization.
        lsVolPointInterp_;
        
        % The number of degrees of freedom of the system.
        degreesOfFreedom_;
        
        vol_x_;
        volDU_;
        volF_
        surfaceF_
        
        input_volU_;
        
        % Projection matrix field: I - nn
        I_nn_;
        
        % For debugging purpose.
        plotErroGradU_ = false;
        isToCheckAgainstAnalyticSolution_ = false;
        
        % Stored for performance reason.
        tensorialIndex_;
    end
    
    methods
        function obj = ...
                NonLinearBlockCoupled(mechanicalModel, problem)
            %NonLinearBlockCoupled Construct an instance of this class
            %   Detailed explanation goes here
            
            obj@SolidModelSolver(mechanicalModel, problem)
            
            obj.lsVolPointInterp_ = ...
                LeastSquareVolPointInterpolation(obj.fvMesh());
            obj.setDegreesOfFreedom();
            obj.preComputeTensorIndices();
            
            obj.vol_x_ = obj.problem.volX();
            obj.volDU_ = obj.volU();
            volI = obj.volGradU.identity();
            obj.volF_ = volI;
            surfaceI = obj.surfaceGradU.identity();
            obj.surfaceF_ = surfaceI;
            
            obj.input_volU_ = obj.volU();
            
            % Projection matrix field: I - nn
            nf = obj.fvMesh().nf();
            obj.I_nn_ = surfaceI - (nf & nf);
        end
        
        %% BEG - Data member access
        function volDU = volDU(obj)
            
            volDU = obj.volDU_;
        end
        
        function volF = volF(obj)
            
            volF = obj.volF_;
        end
        
        function surfaceF = surfaceF(obj)
            
            surfaceF = obj.surfaceF_;
        end
        
        function vol_x_ = vol_x(obj)
           
            vol_x_ = obj.vol_x_;
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
            
            surfacePiola = obj.surfaceF_ * obj.surfaceSigma_;
        end
        
        %function setSurfaceGradU(obj, newValue)
        %    
        %    obj.surfaceGradU_ = newValue;
        %end
        %% END
        
        %% BEG - General member functions.
        function obj = setDegreesOfFreedom(obj)
            
            mesh = obj.fvMesh().spaceMesh();
            obj.degreesOfFreedom_ = 3*(...
                      mesh.numberOfElements ...
                    + mesh.numberOfBElements ...
                    - numberOfFacesInEmptyBoundaries(mesh)...
                    );
        end
        
        function obj = preComputeTensorIndices(obj)
            
            obj.tensorialIndex_ = zeros(obj.degreesOfFreedom_,3);
            
            for i = 1:obj.degreesOfFreedom_
                
                obj.tensorialIndex_(i,:) = getTensorialIndex(i);
            end
        end
        
        % Solve the momentum equation.
        function evolve(obj)
            
            global RUN_OPTIMISATIONS;
            
            fvMesh = obj.fvMesh();
            Sf = obj.fvMesh().Sf();
            mesh = obj.fvMesh().spaceMesh();
            
            % Store previews iteration. Needed for residual computation.
            volU_prevIter = obj.volU();
            
            nCorrectors = obj.config().nCorrectors;
            
            oscillationMonitor = OscillationMonitor(3, 1e-2);
            
            % Calculate a sequence of volDUs (\delta u), the sum will be 
            % the difference x(t_i) - y(t_{i-1}), where y and x are the
            % position at the times t_{i-1} and t_i respectively.
            %
            for iCorrector = 1:nCorrectors
            
                CN = fvMesh.CN();
                gDiff = fvMesh.SfNorm() / CN.norm(); % TODO: add as a mesh prop.

                obj.volDU_ = 0*obj.volDU_;

                % Get second piola stress tensor.
                gradU = obj.surfaceGradU();
                obj.surfaceSigma_ = ...
                    obj.mechanicalModel.computeSigmaField(gradU);
                Sigma = obj.surfaceSigma_;

                nf = obj.fvMesh().nf();
                F = obj.surfaceF_;
                v = F * Sigma * nf;

                % T_d's calculation
                T = {...
                    obj.mechanicalModel.computeSurfaceTd(F,1), ...
                    obj.mechanicalModel.computeSurfaceTd(F,2), ...
                    obj.mechanicalModel.computeSurfaceTd(F,3),...
                };

                % g_d's calculation
                g = {...
                    F.colAsVector(1), ...
                    F.colAsVector(2), ...
                    F.colAsVector(3), ...
                    };

                H = (nf^v)*F.identity()...
                    + (nf^g{1})*T{1}...
                    + (nf^g{2})*T{2}...
                    + (nf^g{3})*T{3};

                % h_d's calculation
                h = {...
                    obj.I_nn_ * g{1}, ...
                    obj.I_nn_ * g{2}, ...
                    obj.I_nn_ * g{3}, ...
                    };
                    
                % Compute first Piola.
                Piola = F * Sigma;
                PSf = Piola * Sf;
                iPSf = PSf.internalField();
                gDiffH = gDiff * H;
                iGDiffH = gDiffH.internalField();
                
                A = zeros(obj.degreesOfFreedom_,obj.degreesOfFreedom_);
                B = zeros(obj.degreesOfFreedom_,1);

                %% BEG - Only for normal derivative and -P(gradU_o)*Sf
                %
                % Beg - Contribution from internal faces.
                %
                for iFace = 1:mesh.numberOfInteriorFaces

                    face = mesh.faces(iFace);

                    coeff = iGDiffH.range(iFace);

                    % beg - Add contribution to matrix diagonal. 
                    i = obj.tensorialIndex_(face.iOwner, :);
                    
                    A(i, i) = A(i, i) - coeff;

                    j = obj.tensorialIndex_(face.iNeighbour, :);
                    A(j, j) = A(j, j) - coeff;
                    % end

                    % Add contribution to matrix off-diagonal.
                    A(i, j) = A(i, j) + coeff;
                    A(j, i) = A(j, i) + coeff;

                    % Computation of the face traction using last
                    % displacement gradient := -P(gradU_o)*Sf
                    B(i) = B(i) - iPSf.range(iFace);
                    B(j) = B(j) + iPSf.range(iFace);
                end
                % End - Contribution from internal faces.


                % Beg - Contribution from boundary faces.
                % 
                for iBoundary=1:mesh.numberOfBoundaries

                    boundary = mesh.boundaries(iBoundary);

                    if strcmp(boundary.type,'empty')
                        continue;
                    end

                    startFace = boundary.startFace;
                    bGDiffH = gDiffH.boundaryField(iBoundary).field();
                    bPSf = PSf.boundaryField(iBoundary).field();

                    for iFace = boundary.startFace:boundary.endFace

                        face = mesh.faces(iFace);
                        iFaceLocal = iFace - startFace + 1;

                        coeff = bGDiffH.range(iFaceLocal);

                        % Add contribution to matrix diagonal.
                        i = obj.tensorialIndex_(face.iOwner, :);
                        A(i, i) = A(i, i) - coeff;

                        % Add contribution to matrix off-diagonal.
                        j = obj.tensorialIndex_(face.iBIndex, :);
                        A(i, j) = A(i, j) + coeff;

                        % Computation of the face traction using last
                        % displacement gradient := -P(gradU_o)*Sf
                        B(i) = B(i) - bPSf.range(iFaceLocal);
                    end
                end
                % End - Contribution from boundary faces.
                %
                %% END - Only for normal derivative and -P(gradU_o)*Sf.


                %% BEG - Tangential derivative
                %
                invLsMatrices = obj.lsVolPointInterp_.invLsMatrices();
                weights = obj.lsVolPointInterp_.weights();
                origins = obj.lsVolPointInterp_.origins();
                mPlaneTransf = ...
                    obj.lsVolPointInterp_.mirrorPlaneTransformations();
            
                if RUN_OPTIMISATIONS
                    % Get internal (raw) data to prevent the slow notation
                    % T{i}.internalField(j) by using T{i}(:,:,j) instead.
                    iT = {...
                        T{1}.internalField().data, ...
                        T{2}.internalField().data, ...
                        T{3}.internalField().data,...
                    };
                
                    % Get internal (raw) data to prevent the slow notation
                    % h{i}.internalField(j) by using h{i}(:,j) instead.
                    ih = {...
                        h{1}.internalField().data, ...
                        h{2}.internalField().data, ...
                        h{3}.internalField().data,...
                    };
                
                    % Analogous as above
                    iI_nn = obj.I_nn_.internalField().data;
                    iv = v.internalField().data;
                end

                % For each cell, it computes the tangential derivative.
                for iElement = 1:mesh.numberOfElements

                    iFaces = mesh.elements(iElement).iFaces;
                    l = obj.tensorialIndex_(iElement, :);

                    for i = 1:length(iFaces)

                        iFace = iFaces(i);
                        face = mesh.faces(iFace);

                        isBoundaryFace = face.patchIndex > 0;

                        if isBoundaryFace

                            boundary = mesh.boundaries(face.patchIndex);

                            if strcmp(boundary.type,'empty')
                                continue;
                            end

                            startFace = boundary.startFace();
                            iFaceLocal = iFace - startFace + 1;
                        end

                        signal = 1;

                        if iElement ~= face.iOwner
                            signal = -1;
                        end

                        if isBoundaryFace

                            pIndex = face.patchIndex;

                            iFl = iFaceLocal;
                            T_d = {...
                                T{1}.boundaryField(pIndex).field(iFl),...
                                T{2}.boundaryField(pIndex).field(iFl),...
                                T{3}.boundaryField(pIndex).field(iFl)...
                                };

                            I_nn = obj.I_nn_.boundaryField(pIndex).field(iFl);
                            bv = v.boundaryField(pIndex).field(iFl);
                            v_proj = I_nn * bv;

                            h_d = {...
                                h{1}.boundaryField(pIndex).field(iFl),...
                                h{2}.boundaryField(pIndex).field(iFl),...
                                h{3}.boundaryField(pIndex).field(iFl),...
                                };
                        else
                            for k = 0:2
                                idx = 3-k;
                                    
                                if RUN_OPTIMISATIONS
                                    T_d{idx} = ...
                                        signal*iT{idx}(:,:,iFace);
                                    
                                    h_d{idx} = ih{idx}(:,iFace);
                                else
                                    T_d{idx} = ...
                                        signal*T{idx}.internalField(iFace);
                                    
                                    h_d{idx} = h{idx}.internalField(iFace);
                                end
                            end

                            if RUN_OPTIMISATIONS
                                I_nn = iI_nn(:,:,iFace);
                                v_ = iv(:,iFace);
                            else
                                I_nn = obj.I_nn_.internalField(iFace);
                                v_ = v.internalField(iFace);
                            end
                            
                            v_proj = I_nn*(signal*v_);
                        end

                        for edge=face.edges

                            m = edge.biNormal;
                            m_dot_v = (m.' * v_proj) * eye(3,3);
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
                                        NonLinearBlockCoupled. ...
                                        tangentialDerivativeCoeff...
                                        (...
                                            L, w, M, dr, m_dot_v, m,...
                                            h_d,...
                                            T_d...
                                        );

                                    normNf = ...
                                        norm(mPlaneTransf(iNode).normal);

                                    if normNf > eps

                                        U = mPlaneTransf(iNode).tensor;
                                        coeff = coeff + U*coeff;
                                    end

                                    iPElement = iPointElements(k);

                                    if iPElement == iElement
                                        % Add contribution to matrix
                                        % diagonal 
                                        A(l, l) = A(l, l) + coeff;
                                    else
                                        % Add contribution to matrix
                                        % off-diagonal 
                                        c = obj.tensorialIndex_(iPElement,:);
                                        A(l, c) = A(l, c) + coeff;
                                    end

                                    neiID = neiID + 1;
                                end % For each point-cell neighbour
                            end % For each node
                        end % For each edge
                    end % For each cell face
                end % For each cell, it computes the tangential derivative.
                %% END


                %% BEG - Contribution from boundary condition
                %
                for iBoundary=1:mesh.numberOfBoundaries

                    boundary = mesh.boundaries(iBoundary);
                    boundaryU = obj.input_volU_.boundaryField(iBoundary);

                    % TODO: <<==== HARD - CODED!!!!
                    % I still don't know to to handle it elegantly.
                    if SolidModelSolver. ...
                            isTractionBoundaryCondition(boundaryU)

                        % Attention: Remember that the boundary equation
                        % is multiplied by (-1).

                        boundaryU = boundaryU.update();

                        externalTraction = ...
                            boundaryU.externalTraction().field();

                        startFace = boundary.startFace;
                        bGDiffH = ...
                            gDiffH.boundaryField(iBoundary).field();
                        bPSf = PSf.boundaryField(iBoundary).field();

                        t = 1;
                        for face = mesh.faces(boundary.iFaces)

                            % BEG - Normal derivative computation
                            %
                            iFaceLocal = face.index - startFace + 1;
                            coeff = bGDiffH.range(iFaceLocal);

                            % Add contribution to matrix diagonal 
                            l = obj.tensorialIndex_(face.iBIndex,:);
                            A(l, l) = A(l, l) - coeff;

                            % Add contribution to matrix off-diagonal 
                            c = obj.tensorialIndex_(face.iOwner,:);
                            A(l, c) = A(l, c) + coeff;
                            % end - Normal derivative computation

                            % Computation of the face traction using last
                            % displacement gradient := -P(gradU_o)*Sf
                            B(l) = B(l) + bPSf.range(iFaceLocal);

                            traction = externalTraction.data(:,t);
                            t = t + 1;
                            sf = norm(face.Sf);
                            B(l) = B(l) - traction*sf;
                            % END - Normal derivative computation


                            % BEG - Tangential derivative computation
                            %
                            for k = 0:2
                                
                                idx = 3-k;
                                
                                T_d{idx} = ...
                                    T{idx}.boundaryField(iBoundary). ...
                                    field(iFaceLocal);

                                h_d{idx} = ...
                                    h{idx}.boundaryField(iBoundary). ...
                                    field(iFaceLocal);
                            end

                            I_nn = obj.I_nn_.boundaryField(iBoundary);
                            bv = I_nn * v.boundaryField(iBoundary);
                            v_proj = bv.field(iFaceLocal);

                            for edge=face.edges

                                m = edge.biNormal;
                                m_dot_v = (m.'*v_proj)*eye(3,3);
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
                                            NonLinearBlockCoupled. ...
                                            tangentialDerivativeCoeff...
                                            (...
                                                L, w, M, dr, m_dot_v, m,...
                                                h_d, T_d...
                                            );

                                        Nf = mPlaneTransf(iNode).normal;
                                        if norm(Nf) > eps

                                            U = mPlaneTransf(iNode).tensor;
                                            coeff = coeff + U*coeff;
                                        end

                                        iPElmt = iPointElements(k);

                                        if iPElmt == face.iBIndex
                                            % Add contribution to matrix
                                            % diagonal 
                                            A(l, l) = A(l, l) - coeff;
                                        else
                                            % Add contribution to matrix
                                            % off-diagonal 
                                            c = obj.tensorialIndex_(iPElmt,:);
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
                    elseif SolidModelSolver. ...
                            isDisplacementBoundaryCondition(boundaryU)

                        avgTensor = averageDiagonalTensor(A);
                        diagSign = sum(diag(avgTensor));
                        diagSign = diagSign/norm(diagSign);
                        mag = norm(reshape(avgTensor,1,9));
                        scalarFactor = diagSign*(1/sqrt(3))*mag;
                        Z = eye(3,3)*scalarFactor;
                        t = 1;
                        boundaryU = boundaryU.update();
                        U_b = boundaryU.field().data;
                        U_o = ...
                            obj.volU.boundaryField(iBoundary).field().data;

                        for face = mesh.faces(boundary.iFaces)

                            U = U_b(:,t) - U_o(:,t);
                            t = t + 1;
                            l = obj.tensorialIndex_(face.iBIndex,:);
                            A(l,l) = A(l,l) + Z;
                            B(l) = B(l) + scalarFactor*U;
                        end

                    elseif SolidModelSolver. ...
                        isSymmetryBoundaryCondition(boundaryU)

                        I_ = eye(3);

                        for face = mesh.faces(boundary.iFaces)

                            n = face.nf;
                            nn = n * n.';
                            gDiff = norm(face.Sf) / norm(face.CN);
                            T1 = gDiff*(nn - I_);
                            coeff_b = -nn + T1;

                            l = obj.tensorialIndex_(face.iBIndex,:);
                            A(l,l) = A(l,l) + coeff_b;

                            j = obj.tensorialIndex_(face.iOwner,:);
                            A(l, j) = A(l, j) - T1;
                        end

                    elseif SolidModelSolver. ...
                            isEmptyBoundaryCondition(boundaryU)
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
                %% END - Contribution from boundary condition

                % Add body force
                if obj.problem.hasBodyForce

                    bodyForce = obj.problem.getBodyForce();
                    [r, c] = size(bodyForce.data);
                    dofBodyForce = r*c; 
                    B(1:dofBodyForce) = B(1:dofBodyForce) ...
                        - reshape(bodyForce.data, [dofBodyForce 1]);
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
                        
                        l = obj.tensorialIndex_(iElem,:);
                        dV = obj.mesh.elements(iElem).volume;
                        A(l,l) = A(l,l) - LHS(l,l)*dV;
                        B(l) = B(l) - RHS(l)*dV;
                    end
                end

                obj.volDU_ = ...
                    solveBlockCoupledLinearSystem(obj.volDU_, A, B);
                
                % Temporary
                if obj.plotErroGradU_
                    obj.errorGradDUinternalFaces();
                end
            
                
                %% Beg - Update fields.
                obj.vol_x_ = obj.vol_x_ + obj.volDU_;
                obj.setVolU(obj.vol_x_ - obj.problem.volX());
                
                obj.calculateVolGradU();
                obj.calculateSurfaceGradU();
                
                obj.volF_ = obj.volGradU.identity() + obj.volGradU();
                surI = obj.surfaceGradU().identity();
                obj.surfaceF_ = surI + obj.surfaceGradU();
                
                if obj.isToCheckAgainstAnalyticSolution_
                    
                    obj.showAnalyticResiduals();
                end
                %% End - Update fields.
                
                obj.problem.statistics();
                
                %obj.problem.meshDrawer().run(obj.volU(), obj.volGradU());
                %plotVolVectorField(obj.volDU_);
                
                [hasConverged, residual] = ...
                computeResiduals...
                (...
                    iCorrector,...
                    obj.volU(),...    % Current iter.
                    volU_prevIter,... % Previews iter.
                    obj.volU_o(),...   % Old timestep.
                    1,... %TODO
                    obj.config()...
                );

                if hasConverged %&& iCorrector > 1
                    
                    iCorrector = -1;
                    break;
                end
                
                if isnan(residual) || norm(residual) > 1e+2
                    
                    error('NaN detected.')
                end
                
                if norm(residual) > 1e+2
                    
                    error('Solver diverged')
                end
                
                if false && oscillationMonitor.stop(residual)
                    break
                end

                volU_prevIter = obj.volU();
            end
            
            if iCorrector > 1 && iCorrector == nCorrectors
               disp('.____________________________.')
               disp('|  SOLVER DID NOT CONVERGE!  |')
               disp('\___________________________/')
            end

            obj.setVolU_oo(obj.volU_o());
            obj.setVolU_o(obj.volU());
        end
        
        function calculateGradU(obj,fieldName,callback,params)
            
            nCorrections = tryGetOrDefault(params, 'nCorrections',3);
            tolerance = tryGetOrDefault(params, 'tolerance',1e-7);
            converged = false;
            
            for iCorrection = 1:nCorrections
                
                callback();

                if obj.analyticResidual(fieldName) < tolerance
                    converged = true;
                    break;
                end
            end

            if converged
                status = '';
            else
                status = 'did not';
            end

            fprintf(...
                   ['Analytical residual for ' fieldName ' did '...
                   status 'converged at iteration %d\n'], iCorrection...
               );
        end
        
        function calculateVolGradU(obj)
            
             callback = @() obj.correctVolGradU(false);
             
             if ~obj.isToCheckAgainstAnalyticSolution_
                 callback();
                 return
             end
             
             params = struct();
             
             obj.calculateGradU('volGradU', callback, params);
        end
        
        function calculateSurfaceGradU(obj)

            callback = @() obj.correctSurfaceGradU();
            
            if ~obj.isToCheckAgainstAnalyticSolution_
                 callback();
                 return
             end
            
            params = struct();
            
            obj.calculateGradU('surfaceGradU', callback, params);
        end
        
        function showAnalyticResiduals(obj)
            
            obj.analyticResidual('vol_x');
            obj.analyticResidual('volU');
            obj.analyticResidual('surfaceF');
            obj.analyticResidual('volF');
        end
        %% END - General member functions.
        
        function errorGradDUinternalFaces(obj)
            
            mesh = obj.fvMesh().spaceMesh();

            if obj.fvMesh.timeMesh.isFirstTimeStep()
                t_o = obj.fvMesh.timeMesh.startTime();
            else
                t_o = obj.fvMesh.timeMesh.valueOld();
            end
 
            t_ = obj.fvMesh.timeMesh.value();
            F_ = obj.problem.analyticalSurfaceF(t_);
            F_o = obj.problem.analyticalSurfaceF(t_o);
            gradDU = (F_ - F_o)*F_o.inv();
            
%             U_ = obj.problem.analyticalVolU(t_);
%             U_o = obj.problem.analyticalVolU(t_o);
%             dU = U_ - U_o;
            
            iGrad = gradDU.internalField();
            
            residual = 0;
            skippeds = 0;
            maxError = -inf;
            minError = inf;
            nInteriorFaces = mesh.numberOfInteriorFaces;
            plotError = true;
            if plotError
                N = int32(sqrt(double(nInteriorFaces)));
                dx = 1 / double(N);
                dy = 1 / double(N);
                [xq,yq] = meshgrid(...
                        0:dx:1,...
                        0:dy:1 ...
                    );
                X = zeros(1,nInteriorFaces);
                Y = zeros(1,nInteriorFaces);
                W = zeros(1,nInteriorFaces);
            end
            for iFace = 1:nInteriorFaces
 
                face = mesh.faces(iFace);
                dU_P = obj.volDU_.internalField(face.iOwner);
                dU_N = obj.volDU_.internalField(face.iNeighbour);
                nTerm =  (dU_N - dU_P)/(face.Sf.'*face.CN);
                aTerm = iGrad(iFace)*face.nf;
                
                aTermNorm = norm(aTerm);
                
                if plotError
                    x = mesh.faces(iFace).centroid(1);
                    y = mesh.faces(iFace).centroid(2);
 
                    X(iFace) = x;
                    Y(iFace) = y;
                end
                
                if aTermNorm < eps
                    skippeds = skippeds + 1;
                    
                    if plotError
                        W(iFace) = 0;
                    end
                    
                    continue;
                end
                
                error_ = norm(nTerm - aTerm) / aTermNorm;
                
                if error_ > maxError
                    maxError = error_;
                end
                
                if error_ < minError
                    minError = error_;
                end
                    
                residual = residual + error_;
                
                if maxError > 0.2 && ~plotError
                    error('Exceeded tolerance: %0.2g%%', maxError*100);
                end
                
                if plotError
                    W(iFace) = error_*100;
                    %W(iFace) = norm(nTerm(3));
                end
            end
            
            disp('=============================')
            disp('Error for gradDU*Sf:');
            meanErrorNormilized = ...
                residual/(double(nInteriorFaces)-skippeds)*100;
            sError = sprintf('%0.2g%%', meanErrorNormilized);
            fprintf('Mean error U = %s\n', sError);
            
            maxErrorNormalized = maxError*100;
            sMaxError = sprintf('%0.2g%%', maxErrorNormalized);
            fprintf('Max error U = %s\n', sMaxError);
            
            sMinError = sprintf('%0.2g%%', minError*100);
            fprintf('Min error U = %s\n', sMinError);
            disp('=============================')
            
            if plotError
                vq = griddata(X,Y,W,xq,yq,'nearest');
                %[X_,Y_] = meshgrid(X,Y);
                %maxL = max(W); minL = min(W);
                %levels = [minL:(maxL-minL)/6:maxL];
                [C, h] = contour(xq, yq, vq, 'ShowText','on');
                colorbar;
                %set(h,'LineColor','black')
                %axis equal;
                %axis square;
            end
        end
        
        function [LHS, RHS] = d2dt2(obj, density, U_o, U_oo, problem)
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
                deltaT*U_o.internalField() ...
                - deltaT*U_oo.internalField() ...
                );
            RHS = reshape(right.data, [dof 1]);
            %% END
        end
    end
    
    methods(Static)
        
        function value = M_xx_dF_x_Sf_generic(N, F_o, Sf, dF)

            M = F_o * N * F_o.';
            value = M.doubleDot(dF)*Sf;
        end
        
        function coeff = tangentialDerivativeCoeff(L,w,M,dr,m_dot_v,m,h,T)
            
            weight = w + (M.' * dr);
            
            coeff = ...
                0.5*L*weight*...
                (...
                    m_dot_v ...
                    + (m.'*h{1})*T{1} + (m.'*h{2})*T{2} + (m.'*h{3})*T{3}...
                );
        end
    end
end

