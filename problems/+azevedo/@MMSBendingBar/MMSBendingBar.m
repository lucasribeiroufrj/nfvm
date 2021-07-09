classdef MMSBendingBar < MMSProblem & handle
    %MMSBendingBar See [7, Kamojjala - 2013].
    %
    % TODO:
    %   - Traction-based boundary condition.
    %
    % Nomenclature.: 
    % - The letters U stand for displacement vector
    % respectively.
    % - The U_o, U_oo stands for U at old, and old old time respectively.
    %
    % Ref.:
    % [Kamojjala - 2013] - Verification tests in solid mechanics.
    
    properties
        
        % Max bending angle.
        maxAngle_;
        
        % To draw to mesh.
        meshDrawer_;
        
        % Boundary conditions for displacement U.
        boundaryConditions_;
        
        % Amplitude function.
        B_;
        
        % dBeta/dt.
        dBdt_;
        
        % d2Beta/dt2.
        d2Bdt2_
        
        % The hight of the bar.
        barHeight_;
        
        % Bool to indicate whether to use only Dirichlet boundary cond. or
        % not.
        useTraction_;
        
        % Peak amplitude
        A_;
    end
    
    methods
        
        function obj = MMSBendingBar(casePath, params)
                        
            nTimeSteps = ...
                tryGetOrDefault(params, 'numberOfTimeSteps', 10);
            
            endTime = ...
                tryGetOrDefault(params, 'endTime', 0.5);
            
            period = ...
                tryGetOrDefault(params, 'period', 1.0);
            
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                    'nTimeSteps', nTimeSteps,...
                    'startTime', 0.0,...
                    'endTime', endTime...
                );
            
            obj@MMSProblem(setup, params);
            
            obj.A_ = tryGetOrDefault(params, 'peakAmplitude', pi/2);

            obj.B_ = ...
                tryGetOrDefault...
                (...
                    params,...
                    'beta',...
                    @(time) (obj.A_/2)*(1 - cos(2*pi*time/period)) ...
                );
            
            obj.dBdt_ = ...
                tryGetOrDefault...
                (...
                    params,...
                    'devBeta',...
                    @(time) (obj.A_*pi/period)*(sin(2*pi*time/period))...
                );
            
            obj.d2Bdt2_ = ...
                tryGetOrDefault...
                (...
                    params,...
                    'devBeta',...
                    @(time) (2*obj.A_*(pi/period)^2)...
                        *(cos(2*pi*time/period))...
                );
            
            obj.barHeight_ = tryGetOrDefault(params, 'barHeight', 1.0);
            
            obj.useTraction_ = ...
                tryGetOrDefault(params, 'useTraction', false);
            
            obj.setBoundaryConditions();
            
            obj.setMeshDrawer();
            
            obj.setSolidModel(params);
        end
        
        function obj = setBoundaryConditions(obj)
            
            for iBoundary=1:obj.mesh().numberOfBoundaries
    
                boundary = obj.mesh().boundaries(iBoundary);
                bName = boundary.userName;
                    
                if strcmp(bName,'top') || strcmp(bName,'left') ...
                        || strcmp(bName,'right') || strcmp(bName,'front')...
                        || strcmp(bName,'back')

                    if obj.useTraction_
                        
                        %error('To be implemented.');
                        kind = 'MMSTraction';
                    else
                        kind = 'MMSDisplacement';
                    end
                elseif strcmp(bName,'bottom')

                    % The base is fixed on the ground.
                    kind = 'fixedDisplacement';
                    %kind = 'solidSymmetry';
                else
                    ME = MException('MMSBendingBar:boundaryCondition',...
                        'Invalid boundary name: %s', bName);
                    throw(ME)
                end
                
                obj.boundaryConditions_{iBoundary}.index = iBoundary;
                obj.boundaryConditions_{iBoundary}.problem = obj;
                obj.boundaryConditions_{iBoundary}.kind = kind;
            end
        end
        
        function x = computeVectorField_x(obj, time, X)
            %computeVectorField_x Computation of the position x.
            %   computeVectorField_x(time, iX) returns x as
            %   vectorField at a specific time. 
            %   position: must be a vectorField.
            
            B = obj.B_(time);
            dataX = X.data;
            
            if B < eps
                
                x = vectorField(dataX);
                return
            end
            
            nElements = X.numberOfElements;
            data = zeros(3,nElements);
            H = obj.barHeight_;
            
            for iElement = 1:nElements
                
                X1 = dataX(1,iElement);
                X2 = dataX(2,iElement);
                X3 = dataX(3,iElement);
                    
                R = H / B;
                S = X2 / R;
                x = [...
                        (R + X1)*cos(S) - R;
                        (R + X1)*sin(S);
                        X3...
                    ];
                
                data(:,iElement) = x;
            end
            
            x = vectorField(data);
        end
        
%         function value = displacementAtBoundary(obj, time, iBoundary)
%             
%             mesh = obj.mesh();
%             t = time;
%             boundary = obj.volX.boundaryField(iBoundary);
%             nFaces = boundary.numberOfFaces();
%             interval = boundary.startFace():boundary.endFace();
%             data = zeros(3,nFaces);
%             
%             iFace = 1;
%             for face = mesh.faces(interval)
%             
%                 X = face.centroid;
%                 X1 = X(1);
%                 X2 = X(2);
%                 H = obj.barHeight_;
%                 B = obj.B_(t);
%                 R = H / B;
%                 S = X2 / R;
%                 x = [...
%                         (R + X1)*cos(S) - R;
%                         (R + X1)*sin(S);
%                         X(3)
%                     ];
%                 
%                 % x = F*X
%                 % F = R*U;
%                 % du = x - X == F*X - X == (F-I)*X;
%                 data(:,iFace) = x - X;
%                 
%                 iFace = iFace + 1;
%             end
%             
%             du = vectorField(data);
%             value = patchVectorField(mesh, iBoundary, du);
%         end
        
        function obj = setMeshDrawer(obj)
           
            obj.meshDrawer_ = BoundaryMeshDrawer(obj.mesh());
        end
        
        function U = U(obj)
            %U Displacement field
            
            U = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function U_o = U_o(obj)
            %Uo Displacement field at time "old"
            
            U_o = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function U_oo = U_oo(obj)
            %Uoo Displacement field at time "old-old"
            
            U_oo = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function F = computeTensorFieldF(obj, time, positions)
            %computeTensorFieldF Computation of the deformation F.
            %   computeTensorFieldF(time, numberOfElements) returns F as
            %   tensorField at a specific time.
            %   position: must be a vectorField.
            
            nElements = positions.numberOfElements;
            data = zeros(3,3,nElements);
            I = eye(3,3);

            for iElement = 1:nElements
               
                X = positions.data(:,iElement);
                X1 = X(1);
                X2 = X(2);
                beta = obj.B_(time);
                alpha =  beta*X2 / obj.barHeight_;
                lambda = beta*X1 / obj.barHeight_ + 1;
                cosine = cos(alpha);
                sine = sin(alpha);
                R = [ ...
                        cosine, -sine,  0;...
                        sine,   cosine 0;...
                        0       0       1 ...
                    ];
                U = I;
                U(2,2) = lambda;
                F_ = R*U;
                data(:,:,iElement) = F_;
            end
            
            F = tensorField(data);
        end
        
        function [iF, bF] = analyticalFParts(obj, X, time)
            %analyticalFParts return internal and boundary fields of F.
            %   analyticalFParts(X, time) returns iF and bF at a specific
            %   time as tensorFields.
            %   X: must be a volVectorField or a surfaceVectorField.
            
            mesh = obj.mesh;
            
            % Use X's internalField to compute F's internalField.
            iX = X.internalField();
            iF = obj.computeTensorFieldF(time, iX);

            nBoundaries = mesh.numberOfBoundaries;
            bF = bTensorField(nBoundaries);
            
            for iBoundary = 1:nBoundaries
                
                % Get an X's patch.
                pX = X.boundaryField(iBoundary);
                
                % Use the tensor field from pX to create a patch field.
                field = obj.computeTensorFieldF(time, pX.field());
                pF = patchTensorField(mesh, iBoundary, field);
                
                % Add the calculated patch to F's boundary field.
                bF = bF.setPatch(iBoundary, pF);
            end
        end
        
        function F = analyticalVolF(obj, time)
            %analyticalVolF analytical deformation tensor F.
            %   analyticalVolF(time) returns F at a specific time as a
            %   volTensorField.
            
            % "Mount" the final F field using the parts: 
            X = obj.volX();
            [internalF, boundaryF] = obj.analyticalFParts(X, time);
            F = volTensorField(obj.mesh, internalF, boundaryF);
        end
        
        function F = analyticalSurfaceF(obj, time)
            %analyticalSurfaceF analytical deformation tensor F.
            %   analyticalSurfaceF(time) returns F at a specific time as a
            %   surfaceTensorField.

            % "Mount" the final F field using the parts: 
            X = obj.surfaceX();
            [internalF, boundaryF] = obj.analyticalFParts(X, time);
            F = surfaceTensorField(obj.mesh, internalF, boundaryF);
        end
        
        function meshDrawer = meshDrawer(obj)
            
            meshDrawer = obj.meshDrawer_;
        end
        
        function setCamera(obj)
            
            xLim = 1.2;
            
%             if obj.barHeight_ > 0
%                 xLim = 1 + obj.deltaLength_ + 0.1;
%             end
            
            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],...
                'XLim',[-obj.barHeight_*2 xLim],...
                'YLim',[-0.2 obj.barHeight_*2.0],...
                'ZLim',[-0.2 1.2]...
             )
            xlabel('x') 
            ylabel('y')
            zlabel('z')
            
            camorbit(20, 20,'data',[0 1 0]);
        end
        
        function material = createMaterial...
                (...
                    obj,...
                    materialEnum,...
                    varargin...
                )
            
            material = ...
                materials.createMaterial...
                (...
                    materialEnum,...
                    planeStress(false),...
                    varargin{:}...
                );
        end
        
        function bool = hasBodyForce(obj)
           
            bool = true;
        end
        
        function RHS = bodyForce(obj)
            
            force = obj.getBodyForce();
            
            RHS = -force.data.';
        end
        
        function force = getBodyForce(obj)
            
            %force = obj.getBodyForceFromArticle();
            %return;
            
            mesh = obj.mesh();
            t = mesh.timeMesh.value();
            B = obj.B_(t);
            
            % No deformation implies no force!
            if B < eps
                force = vectorField(zeros(3,mesh.numberOfElements));
                return;
            end
            
            errorMsg = 'Only HookeanElastic, NeoHookean and StVenant are supported';
            
            if isa(obj.materialLaw(), 'materialLaws.NeoHookean')
            
                law = materialLaws.List.NeoHookean;

            elseif isa(obj.materialLaw(), 'materialLaws.HookeanElastic')
                
                law = materialLaws.List.HookeanElastic;

            elseif isa(obj.materialLaw(), 'materialLaws.StVenant')
            
                law = materialLaws.List.StVenant;
            else
                error(errorMsg);
            end
            
            X = obj.volX().internalField();
            nElements = X.numberOfElements;
            data = zeros(3,nElements);
            dataX = X.data;
            Density = obj.materialLaw().density();
            mu = obj.materialLaw().material().mu();
            lambda = obj.materialLaw().material().lambda();
            H = obj.barHeight_;
            
            for iElement = 1:nElements
                
                X1 = dataX(1,iElement);
                X2 = dataX(2,iElement);
                L = 1 + B*X1/H; % Capital lambda
                alpha = B*X2/H;
                cosine = cos(alpha);
                sine = sin(alpha);
                dBdt = obj.dBdt_(t);
                d2Bdt2 = obj.d2Bdt2_(t);
                
                a_r = -(2*H^3*(-1 + cosine) + H*X2^2*B^2 + X1*X2^2*B^3)*dBdt^2/(H^2*B^3) ...
                    + (H^3*(-1 + cosine)*B*d2Bdt2)/(H^2*B^3);
                a_t = (2*H*(H*sine-X2*B))*dBdt^2/(H*B^3) ...
                    +(-H^2*sine + H*X2*B + X1*X2*B^2)*B*dBdt^2/(H*B^3);
                
                %% Beg - TEMPORARY CODE
                I = eye(3,3);
                U = I; %stretchTensor (Articles' U)
                U(2,2) = L;
                J = det(U);
                dJdL = 1;
                dUdL = zeros(3,3);
                dUdL(2,2) = 1;
                
                if law == materialLaws.List.NeoHookean
                    
                    cauchyBar = lambda / J *log(J)*I + mu/J*(U*U - I);
                    
                    dcauchyBardL = (lambda*(1 - log(J)*dJdL)*I ...
                        + mu*(2*U*dUdL*J - (U*U - I)*dJdL))/J^2;
                    
                elseif law == materialLaws.List.HookeanElastic
                    
                    F = U;
                    dFdL = dUdL;
                    E = F - I; % == F - I = gradU == 0.5*(gradU + gradU.')
                    
                    % For small deformation/strain:
                    % cauchyBar == SigmaBar (Snd-Piola);
                    cauchyBar = (lambda*trace(E))*I + (2*mu)*E;
                    
                    % For small deformation/strain:
                    % dcauchyBardL == dSigmaBardl;
                    dEdL = dFdL;
                    dcauchyBardL = (lambda*trace(dEdL))*I + (2*mu)*dEdL;
                    
                elseif law == materialLaws.List.StVenant
                    
                    F = U;
                    dFdL = dUdL;
                    dF_TdL = dUdL;
                    F_T = F.';
                    E = 0.5*(F_T*F - I);
                    recJ = 1/J;
                    SigmaBar = (lambda*trace(E))*I + (2*mu)*E;
                    cauchyBar = 1/J*F*SigmaBar*F_T;
                    
                    dEdL = 0.5*(dF_TdL*F + F_T*dFdL);
                    dSigmaBardl = (lambda*trace(dEdL))*I + (2*mu)*dEdL;
                    dcauchyBardL = -1/(J*J)*dJdL*F*SigmaBar*F_T ...
                        + recJ*dFdL*SigmaBar*F_T ...
                        + recJ*F*dSigmaBardl*F_T ...
                        + recJ*F*SigmaBar*dF_TdL;
                else
                    error(errorMsg)
                end
                
                cauchy11Bar = cauchyBar(1,1);
                cauchy22Bar = cauchyBar(2,2);
                cauchy12Bar = cauchyBar(1,2);
                cauchy21Bar = cauchyBar(2,1);
                
                dcauchy11BardL = dcauchyBardL(1,1);
                dcauchy21BardL = dcauchyBardL(2,1);
                
                K = B/(H*L);
                f1Bar = K*(-cauchy22Bar + L*dcauchy11BardL + cauchy11Bar);
                f2Bar = K*( cauchy12Bar + L*dcauchy21BardL + cauchy21Bar);
                
                % Analytical for NeoHookean
%                 f_r = -B*(...
%                         H^2*lambda*(-1 + log(L)) ...
%                         + 2*H*mu*X1*B ...
%                         + mu*X1^2*B^2 ...
%                     )/(H*(H + X1*B)^2);
                
                density = Density/L;
                f_r = f1Bar;
                f_t = f2Bar;
                e_r = [ cosine; sine ];
                e_t = [ -sine; cosine];
                b_r = a_r - f_r/density;
                b_t = a_t - f_t/density;
                b = b_r*e_r + b_t*e_t;
                data(:,iElement) = [b(1); b(2); 0];
            end
            
            force = vectorField(data)*mesh.volume()*Density;
        end
        
        function force = getBodyForceFromArticle(obj)
            % It seems incorrect. When using it, the errors are very high.
            
            mesh = obj.mesh();
            t = mesh.timeMesh.value();
            B = obj.B_(t);
            
            X = obj.volX().internalField();
            nElements = X.numberOfElements;
            dataX = X.data;
            Density = obj.materialLaw().density();
            mu = obj.materialLaw().mu();
            lambda = obj.materialLaw().lambda();
            H = obj.barHeight_;
            A = obj.A_;
            dataB = zeros(3,nElements);
            
            for iElement = 1:nElements
                
                X1 = dataX(1,iElement);
                X2 = dataX(2,iElement);
                L = 1 + B*X1/H; % Capital lambda
                alpha = B*X2/H;
                
                p1 = 128*H^3 - 8*A^2*H*X2^2 - 5*A^3*X1*X2^2 ...
                    + 4*(16*H^3 + A^2*H*X2^2 + A^3*X1*X2^2)*cos(2*pi*t);
                p2 = 4*A^2*(2*H + A*X1)*X2^2*cos(4*pi*t) ...
                    - 4*A^3*X1*X2^2*cos(6*pi*t) + A^3*X1*X2^2*cos(8*pi*t);
                p3 = -128*H^3*cos(A*X2*sin(pi*t)^2/H) ...
                    -32*H^3*cos(2*pi*t - A*X2*sin(pi*t)^2/H);
                p4 = -32*H^3*cos(2*pi*t + A*X2*sin(pi*t)^2/H);
                p5 = (A*(1 + A*X1*(1 - cos(2*pi*t))/(2*H))) ...
                    /(2*H*Density*(2*H + A*X1 - A*X1*cos(2*pi*t))^2);
                p6 = -8*H^2*lambda + 8*A*H*mu*X1 + 3*A^2*mu*X1^2 ...
                    -4*A*mu*X1*(2*H + A*X1)*cos(2*pi*t);
                p7 = A^2*mu*X1^2*cos(4*pi*t) ...
                    + 8*H^2*lambda*log(1 + A*X1*sin(pi*t)^2/H);
                p8 = -12*A*H*X2 - 4*A^2*X1*X2 ...
                    + A*(8*H + 7*A*X1)*X2*cos(2*pi*t) ...
                    + 4*A*(H - A*X1)*X2*cos(4*pi*t);
                p9 = A^2*X1*X2*cos(6*pi*t) ...
                    + 32*H^2*sin(A*X2*sin(pi*t)^2/H);
                t10 = A*X2*sin(pi*t)^2/H;
                p10 = -8*H^2*sin(2*pi*t - t10) ...
                    + 8*H^2*sin(2*pi*t + t10);
                br = (pi^2*csc(pi*t)^4*(p1 + p2 + p3 + p4))/(32*A*H^2) ...
                    + p5*(p6 + p7);
                bt = (pi^2*csc(pi*t)^4*(p8 + p9 + p10))/(8*A*H);
                bx = br*cos(alpha) - bt*sin(alpha);
                by = br*sin(alpha) + bt*cos(alpha);
                dataB(1:2,iElement) = [bx; by];
            end
            
            force = vectorField(dataB)*mesh.volume()*Density;
        end
        
        function force = tractionOnTheFace(obj)
            % TODO
            
            mesh = obj.mesh();
            t = mesh.timeMesh.value();
            
            B = obj.B_(t);
            
            X = obj.volX().internalField();
            nElements = X.numberOfElements;
            data = zeros(3,nElements);
            dataX = X.data;
            Density = obj.materialLaw().density();
            mu = obj.materialLaw().mu();
            lambda = obj.materialLaw().lambda();
            H = obj.barHeight_;
            dataU = repmat(eye(3,3),1,1,nElements);
            dataQ = repmat(eye(3,3),1,1,nElements);
            
            for iElement = 1:nElements
                
                X1 = dataX(1,iElement);
                X2 = dataX(2,iElement);
                L = 1 + B*X1/H; % Capital lambda
                alpha = B*X2/H;
                cosine = cos(alpha);
                sine = sin(alpha);
                
                dataU(2,2,iElement) = L;
                dataQ(1:2,1:2,iElement) = ...
                    [ ...
                        cos(alpha), -sin(alpha);...
                        sin(alpha), cos(alpha) ...
                    ];
            end
            
            
            U = tensorField(dataU);
            I = mesh.volI().internalField();
            gradU = U - I; % Note that: F == U
            cauchy = obj.solidModel(). ...
                mechanicalModel().materialLaw().computeSigmaField(gradU);
            Q = tensorField(dataQ);
            traction = Q*Cauchy*Normal;
            force = traction;
        end
    end
end

