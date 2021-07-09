classdef BlockCoupledWithNewBoundary < SolidModelSolver & handle
    %BlockCoupledWithNewBoundary Summary of this class goes here
    %   This class implements the algorithm described in [Cardiff - 2016],
    %   but it implements another traction boundary discretization.
    %
    % Ref.:
    % [Cardiff - 2016] - A block-coupled Finite Volume methodology for
    % linear elasticity and unstructured meshes - Cardiff.
    % [Clifford - 2011] - Block-Coupled Simulations Using OpenFOAM - 
    % Clifford
    
    properties
        
        % Deformation gradient field
        volF_;
        
        % Inverse of the deformation gradient field
        volFinv_
        
        % Jacobian of the deformation gradient
        volJ_;
        
        % Implicit stiffness surface field
        surfaceImpK_;
        
        % Implicit stiffness surface field
        surfaceImpKTemp_;
        
        % To compute tangential discretization.
        lsVolPointInterp_;
        
        % The number of degrees of freedom of the system.
        degreesOfFreedom_;
        
        boundariesU_;
    end
    
    methods
        function obj = ...
                BlockCoupledWithNewBoundary(mechanicalModel, problem)
            %BlockCoupledWithNewBoundary Construct an instance of this class
            %   Detailed explanation goes here
            
            if ~strcmp...
                    (...
                        class(mechanicalModel.materialLaw()),...
                        'HookeanElastic'...
                    )
                error('Only HookeanElastic is supported');
            end
            
            obj@SolidModelSolver(mechanicalModel, problem)

            obj.lsVolPointInterp_ = ...
                LeastSquareVolPointInterpolation(obj.fvMesh());
            obj.setDegreesOfFreedom();
            
            nEmptyBoundaries = 0;
            nBoundaries = obj.fvMesh.numberOfBoundaries();
            for iBoundary=1:nBoundaries
                
                boundaryU = obj.volU().boundaryField(iBoundary);
                    
                if SolidModelSolver.isEmptyBoundaryCondition(boundaryU)
                    
                    nEmptyBoundaries = nEmptyBoundaries + 1;
                end
            end

            
            for iBoundary=nBoundaries - nEmptyBoundaries:-1:1
                
                boundaryU = obj.volU().boundaryField(iBoundary);
                obj.boundariesU_{iBoundary} = boundaryU;
            end
        end
        
        %% BEG - Data member access
        function surfaceImpK = surfaceImpK(obj)
            % Implicit stiffness surface field
            
            surfaceImpK = obj.surfaceImpK_;
        end
        
        function obj = setDegreesOfFreedom(obj)
            
            mesh = obj.fvMesh().spaceMesh();
            obj.degreesOfFreedom_ = 3*mesh.numberOfElements;
        end
        %% END
        
        % Solve the momentum equation.
        function evolve(obj)
            
            Sf = obj.fvMesh().Sf();
            mesh = obj.fvMesh().spaceMesh();
            mu = obj.mechanicalModel.materialLaw_.mu();
            lambda = obj.mechanicalModel.materialLaw_.lambda();
            
            % Store previews iteration. Needed for residual computation.
            volU_prevIter = obj.volU_;
            
%             t = obj.fvMesh.timeMesh.value();
%             fieldName='volU';
%             camelName = [upper(fieldName(1)) fieldName(2:end)];
%             avolU_ = obj.problem.(['analytical' camelName])(t);
            
            for iCorrector=1:obj.config().nCorrectors
                
                % Update boundaries
                for iBoundary=1:length(obj.boundariesU_)

                    obj.boundariesU_{iBoundary} = ...
                        obj.volU_.boundaryField(iBoundary).update();
                    %obj.boundariesU_{iBoundary}.update();
                end
                
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
                    
                    coeff = BlockCoupledWithNewBoundary. ...
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
                %
                for iBoundary=1:length(obj.boundariesU_)
                    
                    boundaryU = obj.boundariesU_{iBoundary};
                    
                    boundary = mesh.boundaries(iBoundary);
                    
                    if SolidModelSolver. ...
                            isDisplacementBoundaryCondition(boundaryU)
                        
                        t = 1;
                        for iFace = boundary.startFace:boundary.endFace

                            face = mesh.faces(iFace);
                            
                            coeff = BlockCoupledWithNewBoundary. ...
                                normalDerivativeCoeff(face, mu, lambda);

                            % Add contribution to matrix diagonal.
                            i = getTensorialIndex(face.iOwner);
                            A(i, i) = A(i, i) - coeff;

                            U = boundaryU.data(:,t);
                            t = t + 1;
                            B(i) = B(i) - coeff*U;
                        end
                        
                    elseif SolidModelSolver. ...
                            isTractionBoundaryCondition(boundaryU)
                        
                        externalTraction = boundaryU.externalTraction();
                        
                        t = 1;
                        for face = mesh.faces(boundary.iFaces)
                            
                            c = getTensorialIndex(face.iOwner);
                            traction = externalTraction.data(:,t);
                            t = t + 1;
                            sf = norm(face.Sf);
                            B(c) = B(c) - traction*sf;
                        end
                    else
                        error('Invalid boundary patch.');
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
                            
                            if SolidModelSolver. ...
                                isTractionBoundaryCondition(obj.boundariesU_{face.patchIndex})
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
                                        BlockCoupledWithNewBoundary. ...
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
                                    
                                    if iPElement > mesh.numberOfElements
                                        
                                        bFace = mesh.numberOfInteriorFaces...
                                            + iPElement ...
                                            - mesh.numberOfElements;

                                        iBoundary = mesh.faces(bFace)...
                                            .patchIndex;
                                        
                                        boundary = mesh...
                                            .boundaries(iBoundary);

                                        localIndex = bFace ...
                                            - boundary.startFace + 1;
                                        
                                        u = obj.boundariesU_{iBoundary}...
                                            .field(localIndex);
                                        
                                        B(l) = B(l) - coeff*u;
                                        
                                    elseif iPElement == iElement
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
                %% END - Computes equation [(6), Cardiff - 2016].

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
                        
                        l = getTensorialIndex(iElem);
                        dV = obj.mesh.elements(iElem).volume;
                        A(l,l) = A(l,l) - LHS(l,l)*dV;
                        B(l) = B(l) - RHS(l)*dV;
                    end
                end
                
                obj.volU_ = solveBlockCoupledLinearSystem(obj.volU_, A, B);
                
                toc
%                 obj.volGradU_ = obj.fvc_grad(obj.volU_);
                obj.volGradU_ = ...
                    fvc_grad...
                    (...
                        obj.volU_,...
                        Sf,...
                        obj.config(),...
                        obj.volGradU_,...
                        obj.surfaceGradU()...
                    );
                
                obj.correctSurfaceGradU();
                
%                 t = obj.fvMesh.timeMesh.value();
%                 fieldName='volGradU';
%                 camelName = [upper(fieldName(1)) fieldName(2:end)];
%                 obj.volGradU_ = obj.problem.(['analytical' camelName])(t);
%                 
%                 fieldName='surfaceGradU';
%                 camelName = [upper(fieldName(1)) fieldName(2:end)];
%                 obj.surfaceGradU_ = obj.problem.(['analytical' camelName])(t);

                obj.volU_ = obj.volU_...
                    .correctBoundaryConditions...
                    (...
                        obj.volU_,...
                        obj.volGradU_,...
                        obj.surfaceGradU()...
                    );
                
%                 obj.volU_ = volVectorField...
%                     (...
%                         obj.fvMesh(),...
%                         obj.volU_.internalField(),...
%                         avolU_.boundaryField() ...
%                     );

                obj.problem.statistics();
                
                [hasConverged, residual] = ...
                    computeResiduals...
                    (...
                        iCorrector,...
                        obj.volU_,...     % Current iter.
                        volU_prevIter,... % Previews iter.
                        obj.volU_o(),...   % Old timestep.
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
        
        function volGradPhi = fvc_grad(obj, phi)
            
            config = obj.config();
            
            if strcmp(config.gradInterpolationScheme, 'Gauss-Green')
                
                Sf = obj.fvMesh().Sf();
                
                for gradIterComputation=1:config.maxGradIterComputation
                    
                    if gradIterComputation == 1
                        
                        mesh = phi.mesh();
                        facePhi = ...
                            surfaceVectorField...
                            (...
                                mesh,...
                                vectorField(mesh.numberOfInteriorFaces),...
                                phi.boundaryField()...
                            );
                        indices = 1:mesh.numberOfInteriorFaces();
                        owners = [mesh.faces(indices).iOwner];
                        neighs = [mesh.faces(indices).iNeighbour];
                        facePhi.internalField(1:end) = ...
                            (phi.internalField(owners) + phi.internalField(neighs))/2;
                    else
                        phi_f_prime = facePhi.internalField();
            
                        iOwners = mesh.iOwners();
                        iNeighs = mesh.iNeighs();

                        grad_C = tensorField(volGradPhi.internalField(iOwners));
                        grad_F = tensorField(volGradPhi.internalField(iNeighs));

                        rf = mesh.r_Cf.internalField();
                        rC = vectorField(mesh.r_C_.internalField(iOwners));
                        rF = vectorField(mesh.r_C_.internalField(iNeighs));

                        iFaces = 1:mesh.numberOfInteriorFaces;

                        facePhi.internalField(iFaces) = phi_f_prime ...
                                + 0.5*(grad_C + grad_F)*(rf - 0.5*(rC + rF));
                    end
                    
                    volGradPhi = fvc_GaussGreenGrad(facePhi, Sf);
                end
            else
                error('Gradient interpolation scheme is not supported!');
            end
        end
    end
    
    methods(Static)
        
        function coeff = normalDerivativeCoeff(face, mu, lambda)
            
            delta = norm(face.E);
            d = norm(face.CN);
            n = face.nf;
            nn = n * n.';
            I = eye(3);

            %coeff = (2*mu + lambda)*(nn)*(delta / d) + ...
            %    mu*(I - nn)*(delta / d);
            
            coeff = (mu + lambda)*(nn)*(delta / d) + ...
                mu*(I)*(delta / d);
        end
        
        function coeff = ...
                tangentialDerivativeCoeff(L, w, M, dr, mu, lambda, mn,mn_t)
            
            weight = w + (M.' * dr);
            coeff = 0.5*L*weight*(mu*mn + lambda*mn_t);
        end
    end
end

