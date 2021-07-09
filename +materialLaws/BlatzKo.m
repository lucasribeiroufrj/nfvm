classdef BlatzKo ...
        < materialLaws.NonLinearBlockCoupledBaseImpl ...
        & materialLaws.SegregatedBaseImpl ...
        & handle
    %BlatzKo Blatz-Ko's law [Blatz - 1962, Brockman - 1986].
    %   This material law can be used with Segregated and NLBC.
    %
    % Ref.:
    % [Blatz - 1962] - Application of finite elastic theory to the 
    % deformation of rubbery materials - Blatz.
    % [Ogden - 2001] - Nonlinear Elasticity: Theory and Applications
    % - Ogden.
    % [Brockman - 1986] - On the use of the Blatz-Ko constitutive model in 
    % nonlinear finite element analysis - Brockman.
    
    properties
        
    end
    
    methods
        
        function obj = BlatzKo(mesh, material)
            % BlatzKo creates a BlatzKo material law.
            %   mesh : must be fvMesh-compatible.
            %   material : must be an BlatzKo-compatible material.
            
            obj@materialLaws.SegregatedBaseImpl(mesh, material);
        end
        
        %% BEG - materialLaws.SegregatedBaseImpl methods implementation
        %
        function Sigma = computeSigmaTensorField(obj, gradU)
            % Computes the second piola stress for a given displacement
            %   gradient field (as tensorField). The equation is described 
            %   in [Brockman - 1986].
            %
            %   computeSigmaTensorField(gradU) returns a tensorField. The
            %       argument gradU must be tensorField too.
            
            I = gradU.identity();
            F = I + gradU;
            C = F.'*F;
            
            % Article's H
            R = materialLaws.BlatzKo.R(C);
            
            mu = obj.material().mu();
            nu = 0.3;
            beta = nu / (1 - 2*nu);
            J3 = scalarField.sqrt(C.det());
            e = -(beta + 2);
            
            Sigma = mu*(I - (J3^e)*R);
        end
        %% END - materialLaws.SegregatedBaseImpl methods implementation
        
        
        %% BEG - NonLinearBlockCoupledBaseImpl methods implementation
        %
        function Td = computeTdTensorField(obj, F, n, t)
            % computeTdTensorField returns Td = FC : nt
            %                                    (2,3)
            % where C is elasticity tensor, n is the normal face and t is
            % any vector. The inputs must be tensor, vector and vector
            % fields respectively.
            
            C = F.'*F;
            R = materialLaws.BlatzKo.R(C); % Article's H
            mu = obj.material().mu();
            nu = 0.3; % TODO
            beta = nu / (1 - 2*nu);
            J3 = C.det().sqrt_();
            e1 = -(beta + 4);
            e2 = (beta + 2);
            a = mu*e2*J3^e1;
            b = 2*mu*J3^(-e2);
            Qbar = materialLaws.BlatzKo.Qbar(C, n, t);
            
            Td = a*((F*(R*n)) & (R*t)) - b*(F*Qbar);
        end
        %% END - NonLinearBlockCoupledBaseImpl methods implementation
    end
    
    methods(Static)
        
        function R = R(C)
            % R(C) rturns the article's H using.
    
            R = C;
            
            C11 = C(1,1);
            C21 = C(2,1);
            C31 = C(3,1);
            C22 = C(2,2);
            C32 = C(3,2);
            C33 = C(3,3);

            R(1,1) = C22*C33 - C32*C32;
            R(1,2) = C31*C32 - C21*C33;
            R(1,3) = C21*C32 - C22*C31;
            R(2,1) = R(1,2);
            R(2,2) = C11*C33 - C31*C31;
            R(2,3) = C21*C31 - C11*C32;
            R(3,1) = R(1,3);
            R(3,2) = R(2,3);
            R(3,3) = C11*C22 - C21*C21;
        end
        
        function Q = Qbar(C, n, t)
            %Qbat(C,F,n,nBar) = Q : nt (the : at 2,4).
            
            % Allocate space for Q;
            Q = C;
            
            C11 = C(1,1);
            C21 = C(2,1);
            C22 = C(2,2);
            C31 = C(3,1);
            C32 = C(3,2);
            C33 = C(3,3);
            
            nt = n & t;
            
            n1t2 = nt(1,2);
            n1t3 = nt(1,3);
            n2t1 = nt(2,1);
            n2t3 = nt(2,3);
            n3t1 = nt(3,1);
            n3t2 = nt(3,2);
            
            t01 = C33*n1t2;
            t02 = C32*n1t3;
            t03 = C31*n2t3;
            t04 = C33*n2t1;
            t05 = C31*n3t2;
            t06 = C32*n3t1;
            t07 = C22*n1t3;
            t08 = C21*n2t3;
            t09 = C32*n1t2;
            t10 = C21*n3t2;
            t11 = C22*n3t1;
            t12 = C32*n2t1;
            t13 = C11*n2t3;
            t14 = C21*n1t3;
            t15 = C11*n3t2;
            t16 = C31*n1t2;
            t17 = C21*n3t1;
            t18 = C31*n2t1;
            
            %Q11 = 0;
            Q12 = t01 - t02 + t03 - t04 - t05 + t06;
            Q13 = t07 - t08 - t09 + t10 - t11 + t12;
            %Q21 = -Q12;
            %Q22 = 0;
            Q23 = t13 - t14 - t15 + t16 + t17 - t18;
            %Q31 = -Q13;
            %Q32 = -Q23;
            %Q33 = 0;
            
            Q(1,1) = 0*Q(1,1);
            Q(1,2) =  Q12;
            Q(1,3) =  Q13;
            Q(2,1) = -Q12;
            Q(2,2) = 0*Q(2,2);
            Q(2,3) =  Q23;
            Q(3,1) = -Q13;
            Q(3,2) = -Q23;
            Q(3,3) = 0*Q(3,3);
        end
        
        function value = leviCivitaTensor()
           
            persistent lcMat
            
            if isempty(lcMat)
                
                % leviCivitaTensor
                lcMat = zeros(3,3,3);
                lcMat([8 12 22]) = 1;
                lcMat([6 16 20]) = -1;
            end
            
            value = lcMat;
        end
        
        function value = greenDeformationTensor()
            % GreenDeformationTensor() return the symbolic C = F^T * F
            
            persistent C;
            
            if isempty(C)
                
                C = sym('C%d%d', [3 3]);
            end
            
            value = C;
        end
        
        function value = symmetricGreenDeformationTensor()
            % GreenDeformationTensor() return the symbolic C = F^T * F
            
            persistent C;
            
            if isempty(C)
                
                C = materialLaws.BlatzKo.greenDeformationTensor();

                % Set as symmetric
                for a = 1:3
                    for b = a:3
                            C(a,b) = C(b,a);
                    end
                end
            end
            
            value = C;
        end
        
        function R = symbolicRComponent(i,j)
            % R(i,j) compute equation [(8), Blatz - 1962]

            lcMat = materialLaws.BlatzKo.leviCivitaTensor();

            C = materialLaws.BlatzKo.symmetricGreenDeformationTensor();

            R = 0;
            for m = 1:3
               for n = 1:3
                   for p = 1:3
                       for q = 1:3
                           ee = lcMat(i,m,p)*lcMat(j,n,q);
                           R = R ...
                               + 1/2*ee*C(m,n)*C(p,q);
                       end 
                   end
               end 
            end
        end
        
        function value = Q(i,j,m,n)
            % Q(i,j) compute equation [(12), Blatz - 1962]

            lcMat = materialLaws.BlatzKo.leviCivitaTensor();
            C = materialLaws.BlatzKo.symmetricGreenDeformationTensor();

            value = 0;
            for p = 1:3
               for q = 1:3
                   ee = lcMat(i,m,p)*lcMat(j,n,q);
                   value = value ...
                       + ee*C(p,q);
               end 
            end
        end
        
        function value = Q2dotsNNBar(i,m)
            % Q2dotsNNBar(i,j) == (Q:nnBar)(i,m) | nBar = F dot n
            
            nM = sym('N%dM%d', [3 3]);
            
            value = 0;
            for j = 1:3
               for n = 1:3
                   value = value ...
                       + materialLaws.BlatzKo.Q(i,j,m,n)*nM(j,n);
               end 
            end
        end
        
        function value = QBar()
            % QBar() == Q:nnBar
            
            value = sym('Q%d%d', [3 3]);
            
            for i = 1:3
                for m = 1:3
                    value(i,m) = materialLaws.BlatzKo.Q2dotsNNBar(i,m);
                end
            end
        end
        
        function R = symbolicR()
            % symbolicR() == R as symbolic matrix
            
            R = sym('Q%d%d', [3 3]);
            
            for i = 1:3
                for m = 1:3
                    R(i,m) = ...
                        materialLaws.BlatzKo.symbolicRComponent(i,m);
                end
            end
        end
        
        function value = HforNonLinearBlockCoupled(i,m)
            % HforNonLinearBlockCoupled(i,m) computes D_{ijmn}*n_j*n_n. The
            % tensor D is the tangent modulus tensor defined in equation 
            % [(12), Blatz - 1962].
            
            %syms mu B J_3
            %a = mu*(B + 2)*J_3^-(B + 4);
            %b = 2*mu*J_3^-(B+2);
            syms a b
            N = sym('n%d', [1 3]); % face normal.
            
            value = 0;
            for j = 1:3
                for n = 1:3
                    H1 = materialLaws.BlatzKo.H(i,j);
                    H2 = materialLaws.BlatzKo.H(m,n);
                    Q = materialLaws.BlatzKo.Q(i,j,m,n);
                    value = value + (a*H1*H2 - b*Q)*N(j)*N(n);
                end
            end
        end
        
        function value = HforNonLinearBlockCoupled2(i,m)
            % HforNonLinearBlockCoupled(i,m) computes D_{ijmn}*n_j*n_n. The
            % tensor D the tangent modulus tensor defined in equation 
            % [(12), Blatz - 1962].
            
            %syms mu B J_3
            %a = mu*(B + 2)*J_3^-(B + 4);
            %b = 2*mu*J_3^-(B+2);
            syms a b
            %N = sym('n%d', [1 3]); % face normal.
            D = sym('D%d%d', [1 3]); % face normal.
            
            value1 = 0;
            value2 = 0;
            for j = 1:3
                for n = 1:3
                    H1 = materialLaws.BlatzKo.H(i,j);
                    H2 = materialLaws.BlatzKo.H(m,n);
                    Q = materialLaws.BlatzKo.Q(i,j,m,n);
                    value1 = value1 + H1*H2;
                    value2 = value2 - Q;
                    D(j,n) = a*(value1) + b*(value2);
                end
            end
            
            value = D;
        end
    end
end

