classdef OgdenStoraker ...
        < materialLaws.NonLinearBlockCoupledBaseImpl ...
        & materialLaws.SegregatedBaseImpl ...
        & handle
    %OgdenStoraker implements the nth order Ogden model [Hackett - 2016]
    %   See all [Schrodt - 2005].
    %
    % Obs.:
    % (1) It is used here Ogden's notation [Hackett - 2016]. Here the 
    %     Nominal is symmmtrical, then N = N^T = first Piola.
    %
    % Ref.:
    % [Hackett - 2016] - Hyperelasticity Primer
    % [Ogden - 1997] - Non-Linear Elastic Deformations.
    % [Schrodt - 2005] - Hyperelastic Description of Polymer Soft Foams 
    % at Finite Deformations
    
    properties
        
    end
    
    methods
        
        function obj = OgdenStoraker(mesh, material)
            % OgdenStoraker creates a OgdenStoraker material law.
            %   mesh : must be fvMesh-compatible.
            %   material : must be an OgdenStoraker-compatible
            %   material.
            
            obj@materialLaws.SegregatedBaseImpl(mesh, material);
        end
        
        function Sigma = computeSigmaField(obj, gradU)
            % computeSigmaField Computes the second-piola stress 
            %   tensor for a given displacement gradient field.
            %   gradU : volTensorField, surfaceTensorField or 
            %   patchTensorField.
            
            % NOTE: this function is needed by both 
            % NonLinearBlockCoupledBaseImpl and SegregatedBaseImpl.
            
            if isa(gradU, 'patchTensorField')
                
                Sigma = obj.computeSigmaPatchField(gradU);
                return;
            end
            
            % Sigma's internal field.
            [iSigma, iImpK] = ...
                obj.computeSigmaTensorField(gradU.internalField());
            
            nBoundaries = gradU.mesh().numberOfBoundaries;
            bSigmaField = bTensorField(nBoundaries);
            bImpKField = bTensorField(nBoundaries);
            mesh = gradU.mesh();
            
            for iBoundary = 1:nBoundaries
                
                % Get gradU's patch.
                pGradU = gradU.boundaryField(iBoundary).field;
                
                % Use the patch to calculate a tensor field.
                [SigmaField, impKField] = ...
                    obj.computeSigmaTensorField(pGradU);
                
                % Use the tensor field to create a patch field.
                pSigma = patchTensorField(mesh, iBoundary, SigmaField);
                pImpK = patchTensorField(mesh, iBoundary, impKField);
                
                % Add the calculated patch to Sigma's boundary field.
                bSigmaField = bSigmaField.setPatch(iBoundary, pSigma);
                bImpKField = bImpKField.setPatch(iBoundary, pImpK);
            end
            
            % "Mount" the final Sigma field using the parts: 
            % internalfield + boundaryField
            if isa(gradU, 'surfaceTensorField')
                
                Sigma = surfaceTensorField(mesh, iSigma, bSigmaField);
                impK = surfaceTensorField(mesh, iImpK, bImpKField).trace();
                
                obj.setSurfaceImpK(impK);
            else
                Sigma = volTensorField(mesh, iSigma, bSigmaField);
                impK = volTensorField(mesh, iImpK, bImpKField).trace();
                
                obj.setVolImpK(impK);
            end
        end
        
        function [Sigma, impK, U1, U2, U3, U11, U22, U33, eigValues] = ...
                computeSigmaTensorField(obj, gradU)
            %computeSigmaTensorField Compute 2nd Piola(1) and impK from
            %   C (= F^T*F right Cauchy-Green). The gradU input must be a
            %   tensorField. 
            %
            %   See [eq. 5.44, Hackett - 2016]
            
            I = gradU.identity();
            F = I + gradU;
            C = F.'*F;
            
            % Capture eigen vector/values and tensor fields.
            [eigVectors, eigValues] = C.eigen();
            
            % Each column of eigVectors(iElement) stores an eigen vector
            % of C.
            U1 = eigVectors.colAsVector(1); % A field with the 1st e.v.'s.
            U2 = eigVectors.colAsVector(2); % A field with the 2nd e.v.'s.
            U3 = eigVectors.colAsVector(3); % A field with the 3st e.v.'s.
            
            nElements = C.numberOfElements();
            
            % 2nd Piola
            Sigmas = zeros(3,3,nElements);
            
            % 2nd derivative of energy function U(eigenValue_i)
            impKs = zeros(3,3,nElements);
            
            % BEG - For each element calculate Td
            %
            % Get raw data to improve performance and use tensorial product
            % to construct a base.
            eigValues_ = eigValues.data;
            U11 = U1 & U1;
            U22 = U2 & U2;
            U33 = U3 & U3;
            U11_ = U11.data;
            U22_ = U22.data;
            U33_ = U33.data;
            
            alphas = obj.material.alphas();
            betas = obj.material.betas();
            mus = obj.material.mus();
            
            for iElem = 1:nElements
                
                % Beg - Calculate 2nd Piola and K
                EIG = eigValues_(:,:,iElem);
                eigV = [ EIG(1,1) EIG(2,2) EIG(3,3) ];
                J = prod(eigV);
                Sigma_ = zeros(3,3);
                K_ = zeros(3,3);
                
                for iEig = 1:3
                
                    lamda = eigV(iEig);
                    diag = 0;
                    
                    for i = 1:3
                        
                        alpha = alphas(i) / 2;
                        beta = betas(i);
                        mu = mus(i);
                        
                        Sigma_(iEig, iEig) = Sigma_(iEig, iEig) ...
                            + (mu/alpha) ...
                            * (lamda^alpha - (J)^(-alpha*beta));
                        
                        K_(iEig,iEig) = K_(iEig,iEig) ...
                            + (mu/alpha) ...
                            * (...
                                 (alpha-1)*lamda^alpha ...
                                 + (1 + alpha*beta)*(J)^(-alpha*beta)...
                              );
                        
                        if iEig == 3
                            
                            diag = diag + mu*beta*(J)^-(alpha*beta + 1);
                        end
                    end
                    
                    Sigma_(iEig, iEig) = Sigma_(iEig, iEig) / lamda;
                    K_(iEig,iEig)  = K_(iEig,iEig)*(1/(2*lamda*lamda));
                end
                
                K_(1,2) = (eigV(3)/2)*diag;
                K_(2,1) = K_(1,2);
                K_(1,3) = (eigV(2)/2)*diag;
                K_(3,1) = K_(1,3);
                K_(2,3) = (eigV(1)/2)*diag;
                K_(3,2) = K_(2,3);
                
                Sigmas(:,:,iElem) = ...
                      Sigma_(1,1)*U11_(:,:,iElem) ...
                    + Sigma_(2,2)*U22_(:,:,iElem) ...
                    + Sigma_(3,3)*U33_(:,:,iElem);
                impKs(:,:,iElem) = K_;
                % End - Calculate 2nd Piola and K
            end
            
            Sigma = tensorField(Sigmas);
            impK = tensorField(impKs);
        end
        
        
        %% BEG - NonLinearBlockCoupledBaseImpl methods implementation
        %
        function Td = computeTdTensorField(obj, F, n, t)
            % computeTdTensorField returns Td = FC : nt
            %                                    (2,3)
            % where C is elasticity tensor, n is the normal face and t is
            % any vector. The inputs must be tensor, vector and vector
            % fields respectively.
            
            C = F.'*F; % Right Cauchy-Green tensor
            
            nElements = C.numberOfElements();

            [Sigma, K, U1, U2, U3, U11, U22, U33, eigValues] = ...
                obj.computeSigmaTensorField(C);
            
            % Use tensorial product to construct a base and get raw data to
            % improve performance.
            U12 = (F*U1) & U2;
            U13 = (F*U1) & U3;
            U21 = (F*U2) & U1;
            U23 = (F*U2) & U3;
            U31 = (F*U3) & U1;
            U32 = (F*U3) & U2;
            U = {...
                U11.data, U12.data, U13.data; ...
                U21.data, U22.data, U23.data; ...
                U31.data, U32.data, U33.data
                };

            % See above
            Td = tensorField(nElements);
            
            % Pre-compute dot products and get raw data to improve
            % performance.
            U1n = U1 ^ n;
            U2n = U2 ^ n;
            U3n = U3 ^ n;
            U1t = U1 ^ t;
            U2t = U2 ^ t;
            U3t = U3 ^ t;
            Udotn = { U1n.data, U2n.data, U3n.data };
            Udott = { U1t.data, U2t.data, U3t.data };
            Sigma_ = Sigma.data;
            
            % BEG - For each element calculate Td
            for iElem = 1:nElements
                
                EIG = eigValues(iElem);
                eigV = [ EIG(1,1) EIG(2,2) EIG(3,3) ];
                K_ = K(iElem);
                
                % First parcel of Td
                X = zeros(3,3);
                for a = 1:3
                    
                    Uadn = Udotn{a}(iElem);
                    
                    for b = 1:3
                        
                        Ubdt = Udott{b}(iElem);
                        Uab = U{a,b}(:,:,iElem);
                        X = X + K_(a,b)*(Uadn*Ubdt)*Uab;
                    end
                end
            
                % Second parcel of Td
                Y = zeros(3,3);
                for a = 1:3
                    for b = 1:3
                        
                        if a == b
                            continue;
                        end
                        
                        la = eigV(a);
                        lb = eigV(b);
                        
                        if abs(la - lb) > eps
                            
                            iSigma = Sigma_(:,:,iElem);
                            term = (iSigma(a,a)-iSigma(b,b)) / (la-lb);
                        else
                            term = 2*(K_(b,b) - K_(a,b));
                        end
                        
                        Ubdn = Udotn{b}(iElem);
                        Ubdt = Udott{b}(iElem);
                        Uadt = Udott{a}(iElem);
                        Uaa = U{a,a}(:,:,iElem);
                        Uab = U{a,b}(:,:,iElem);
                            
                        Y = Y + term*((Ubdn*Ubdt)*Uaa + (Ubdn*Uadt)*Uab);
                    end
                end
                
                Td(iElem) = 4*X + Y;
            end
            % END - For each element calculate Td
        end
        %% END - NonLinearBlockCoupledBaseImpl methods implementation
    end
end

