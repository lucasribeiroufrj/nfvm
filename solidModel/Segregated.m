classdef Segregated < SolidModelSolver & handle
    %Segregated This is the "standard" cell-centred FVM.
    %   For more information see [Cardiff - 2018]. It uses the Total
    % Lagrangian Total Displacement framework.
    %
    % [Cardiff - 2018] Thirty years of the finite volume method for solid
    % mechanics.
    
    properties (Access = private)
        
        % First piola stress tensor
        volPiola_;
        
        % First piola stress tensor
        surfacePiola_;
        
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
    end
    
    methods
        function obj = Segregated(mechanicalModel, problem)
                %Segregated(mechanicalModel, problem, config)
            
            obj@SolidModelSolver(mechanicalModel, problem)
            
            % Force ImpK calculation
            obj.mechanicalModel().computePiolaField(obj.surfaceGradU());
            obj.surfaceImpK_ = obj.mechanicalModel().surfaceImpK();
        end
        
        %% BEG - Data member access
        function volPiola = volPiola(obj)
            %volPiola First piola stress tensor.
            %   volPiola() returns the first piola stress volume tensor.
            
            volPiola = obj.volPiola_;
        end
        
        function surfacePiola = surfacePiola(obj)
            %surfacePiola First piola stress tensor.
            %   surfacePiola() returns the first piola stress surface
            %   tensor.
            
            surfacePiola = obj.surfacePiola_;
        end
        
        function surfaceImpK = surfaceImpK(obj)
            % Implicit stiffness surface field
            
            surfaceImpK = obj.surfaceImpK_;
        end
        %% END
        
        % Solve the momentum equation.
        function evolve(obj)
            
            Sf = obj.fvMesh().Sf();
            
            % Store previews iteration. Needed for residual computation.
            volU_prevIter = obj.volU();
            
            for iCorrector=1:obj.problem().configuration().nCorrectors
                
                fprintf...
                (...
                    'Iteration %d/%d\n',...
                    iCorrector,...
                    obj.problem().configuration().nCorrectors...
                );
                
                obj.updateBoundaryAndFields();
                
                % Obs.: This code will be organised once we implement the
                % sparse matrix class.
                [LHS, RHS] = ...
                    fvm_laplacian... % equals to divImpKgradU
                    (...
                        obj.surfaceImpK(),...
                        obj.surfaceGradU(),...
                        obj.volU()...
                    );
                
                A = LHS;
                B = RHS + fvc_div(obj.surfacePiola(), Sf)...
                        - fvc_laplacian... % equals to divImpKgradU
                          (...
                            obj.surfaceImpK(),...
                            obj.surfaceGradU(),...
                            Sf...
                          );
                      
                if obj.problem().hasBodyForce
                    B = B + obj.problem().bodyForce();
                end
                      
                if obj.problem.includeAcceleration()
                    
                    [LHS, RHS] = ...
                        fvm_d2dt2...
                        (...
                            obj.mechanicalModel().density(),...
                            obj.volU_o(),...
                            obj.volU_oo(),...
                            obj.problem()...
                        );

                    A = A - LHS;
                    B = B - RHS;
                end
                
                % == 0; i.e. all the term are on the left-hand side of the 
                % equation.
                
                
                % Solve the linear sytem and save the solution in U's
                % internal field.
                %
                obj.setInternalDataVolU(...
                    solveLinearSystem...
                    (...
                        obj.volU().internalField().data,...
                        A,...
                        B,...
                        obj.problem().configuration()...
                    ));
                
                obj.correctVolGradU();
                obj.correctSurfaceGradU();

                obj.problem().statistics();
                
                [hasConverged, residual] = ...
                    computeResiduals...
                    (...
                        iCorrector,...
                        obj.volU(),...     % Current iter.
                        volU_prevIter,...  % Previews iter.
                        obj.volU_o(),...   % Old timestep.
                        1,... %TODO
                        obj.problem().configuration()...
                    );

                if hasConverged && iCorrector > 1
                    
                    break;
                end
                
                if isnan(residual)
                    
                    error('NaN detected.')
                end
                
                volU_prevIter = obj.volU();
            end
            
            % Update fields
            obj.setVolU_oo(obj.volU_o());
            obj.setVolU_o(obj.volU());
            obj.surfaceImpK_ = obj.surfaceImpKTemp_;
        end
        
        function divImpKgradU = ...
                tractionBoundaryLaplacianOfimpKGradU...
                (...
                    obj,...
                    externalTraction,...
                    patch...
                )
            %tractionBoundaryLaplacianImpkDisp Computes laplacian(K*gradU)
            %   In other words, it computes div(K*gradU). Prescibed 
            % solidtraction boundary condition uses it.

            idx = patch.index();
            extTraction = externalTraction;
            S = obj.fvMesh().Sf().boundaryField(idx);
            
            P = obj.surfacePiola().boundaryField(idx);
            impK = obj.surfaceImpK_.boundaryField(idx);
            gradU = obj.surfaceGradU().boundaryField(idx);
            internalForce = P*S;
            I = obj.fvMesh().surfaceI();
            I = I.boundaryField(idx);
            F = I + gradU;
            invF = F.inv();
            J = F.det();
            currentS = J*(invF.' * S);
            externalForce = extTraction*currentS.norm();

            divImpKgradU = externalForce - internalForce + impK*(gradU*S);
        end
    end
    
    methods (Access = private)
        
        function updateBoundaryAndFields(obj)

            surfacePiola = ...
                obj.mechanicalModel().computePiolaField(obj.surfaceGradU());
            
            obj.surfaceImpKTemp_ = obj.mechanicalModel().surfaceImpK();

            obj.setSurfacePiola(surfacePiola);
            obj.setVolU(obj.volU().updateDivImpKgradUfluxes(obj));
        end
    end
    
    methods (Access = protected)
        
       function obj = setVolPiola(obj, newValue)
            % setVolPiola Sets the first piola stress.
            %   setVolPiola(newValue) sets the first piola stress volume
            %   tensor using newValue as a volTensorField.
            
            obj.volPiola_ = newValue;
        end

        function obj = setSurfacePiola(obj, newValue)
            % setSurfacePiola Sets the first piola stress.
            %   setSurfacePiola(newValue) sets the first piola stress
            %   surface tensor using newValue as a surfaceTensorField.
            
            obj.surfacePiola_ = newValue;
        end 
    end
end

