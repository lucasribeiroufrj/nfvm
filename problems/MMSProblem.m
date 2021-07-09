classdef MMSProblem < Problem & handle
    %MMSProblem Base class for problems implementing the Method of
    %   Manuctured solutions.
    %
    %   The deformation gradient F is cosidered as the primary unknown.
    
    properties
        
        % A "table" to store the Mean, Max. and Min. Errors.
        displacementErrors_;
    end
    
    methods
        function obj = MMSProblem(setup, params)
            
            if nargin == 1
                params = [];
            end
            
            obj@Problem(setup, params);
            
            obj.displacementErrors_ = zeros(setup.nTimeSteps,4);
        end
    end
    
    methods (Abstract)
        
        F = computeTensorFieldF(obj, time, X)
        %computeTensorFieldF Computation of the deformation F.
        %   computeTensorFieldF(time, X) returns F as
        %   tensorField at a specific time.
        %   X (position) must be a vectorField.
        
        x = computeVectorField_x(obj, time, X)
        %computeVectorField_x Computation of the position x.
        %   computeVectorField_x(time, X) returns x as
        %   vectorField at a specific time. 
        %   X (position) must be a vectorField.
    end
    
    methods
        
        function [iF, bF] = analyticalF_Parts(obj, X, time)
            %analyticalF_Parts return internal and boundary fields of F.
            %   analyticalF_Parts(X, time) returns iF and bF at a specific
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
            [internalF, boundaryF] = obj.analyticalF_Parts(X, time);
            F = volTensorField(obj.mesh, internalF, boundaryF);
        end
        
        function F = analyticalSurfaceF(obj, time)
            %analyticalSurfaceF analytical deformation tensor F.
            %   analyticalSurfaceF(time) returns F at a specific time as a
            %   surfaceTensorField.

            % "Mount" the final F field using the parts: 
            X = obj.surfaceX();
            [internalF, boundaryF] = obj.analyticalF_Parts(X, time);
            F = surfaceTensorField(obj.mesh, internalF, boundaryF);
        end
        
        function [ix, bx] = analytical_x_Parts(obj, X, time)
            %analytical_x_Parts return internal and boundary fields of x.
            %   analytical_x_Parts(X, time) returns ix and bx at a specific
            %   time as tensorFields.
            %   X: must be a volVectorField or a surfaceVectorField.
            
            mesh = obj.mesh;
            
            % Use X's internalField to compute F's internalField.
            iX = X.internalField();
            ix = obj.computeVectorField_x(time, iX);
            
            nBoundaries = mesh.numberOfBoundaries;
            bx = bVectorField(nBoundaries);
            
            for iBoundary = 1:nBoundaries
                
                % Get an X's patch.
                pX = X.boundaryField(iBoundary);
                
                % Use the tensor field from pX to create a patch field.
                field = obj.computeVectorField_x(time, pX.field());
                px = patchVectorField(mesh, iBoundary, field);
                
                % Add the calculated patch to x's boundary field.
                bx = bx.setPatch(iBoundary, px);
            end
        end
        
        function x = analyticalVol_x(obj, time)
            %analyticalVol_x analytical centroid position x.
            %   analyticalVol_x(time) returns x at a specific time as a
            %   volVectorField.
           
            mesh = obj.mesh;
            
            % "Mount" the final F field using the parts: 
            X = obj.volX();
            [internal_x, boundary_x] = obj.analytical_x_Parts(X, time);
            x = volVectorField(mesh, internal_x, boundary_x);
        end
        
        function x = analyticalSurface_x(obj, time)
            %analyticalSurface_x analytical centroid position x.
            %   analyticalSurface_x(time) returns x at a specific time as a
            %   surfaceVectorField.
           
            mesh = obj.mesh;
            
            % "Mount" the final F field using the parts: 
            X = obj.surfaceX();
            [internal_x, boundary_x] = obj.analytical_x_Parts(X, time);
            x = surfaceVectorField(mesh, internal_x, boundary_x);
        end
        
        function U = analyticalVolU(obj, time)
            %analyticalVolU analytical deformation U.
            %   analyticalVolU(time) returns U at a specific time as a
            %   volVectorField.
            
            x = obj.analyticalVol_x(time);
            X = obj.volX();
            U = x + (-1)*X;
        end
        
        function U = analyticalSurfaceU(obj, time)
            %analyticalSurfaceU analytical deformation U.
            %   analyticalSurfaceU(time) returns U at a specific time as a
            %   surfaceVectorField.
            
            x = obj.analyticalSurface_x(time);
            X = obj.surfaceX();
            U = x + (-1)*X;
        end
        
        function gradU = analyticalSurfaceGradU(obj, time)
            %analyticalSurfaceGradU analytical gradient gradU.
            %   analyticalSurfaceGradU(time) returns gradU at a specific
            %   time as a surfaceTensorField.
            
            I = obj.mesh().surfaceI();
            F = obj.analyticalSurfaceF(time);
            gradU = F - I;
        end
        
        function gradU = analyticalVolGradU(obj, time)
            %analyticalVolGradU analytical gradient gradU.
            %   analyticalVolGradU(time) returns gradU at a specific
            %   time as a volTensorField.
            
            I = obj.mesh().volI();
            F = obj.analyticalVolF(time);
            gradU = F - I;
        end
        
        function Sf = analyticalDeformedSf(obj, time)
            %analyticalDeformedSf analytical surface vector Sf.
            %   analyticalDeformedSf(time) returns Sf at a specific
            %   time as a surfaceVectorField.
            
            Sf = obj.deformedSf(obj.analyticalSurfaceF(time));
        end
        
        function Sf = deformedSf(obj, F)
            
            inv_F = F.inv();
            J = F.det();
            mesh = obj.mesh();
            S = ...
                surfaceVectorField...
                (...
                    mesh,...
                    vectorField([ mesh.faces.Sf ])...
                );
            Sf = J*(inv_F.' * S);
        end
        
%         function obj = displacementError(obj, phi, time)
%             
%             mesh = obj.mesh();
%             U = obj.analyticalVolU(time);
%             N = mesh.numberOfElements;
%             
%             residual = 0;
%             skippeds = 0;
%             maxError = -inf;
%             minError = inf;
%             for iElem = 1:N
%                 nU = phi.internalField(iElem);    % numerical U
%                 aU = U.internalField(iElem);      % analytical U
%                 P = norm(aU);
%                 
%                 if P < eps
%                     skippeds = skippeds + 1;
%                     continue;
%                 end
%                 
%                 error = norm(nU - aU) / P;
%                 
%                 if error > maxError
%                     maxError = error;
%                 end
%                 
%                 if error < minError
%                     minError = error;
%                 end
%                     
%                 residual = residual + error;
%             end
%             
%             meanErrorPerc = residual/(N-skippeds)*100;
%             sError = sprintf('%0.2g%%', meanErrorPerc);
%             fprintf('Mean error U = %s\n', sError);
%             
%             maxErrorPerc = maxError*100;
%             sMaxError = sprintf('%0.2g%%', maxErrorPerc);
%             fprintf('Max  error U = %s\n', sMaxError);
%             
%             minErrorPerc = minError*100;
%             sMinError = sprintf('%0.2g%%', minErrorPerc);
%             fprintf('Min  error U = %s\n', sMinError);
%             
%             iTimeStep = mesh.timeMesh.timeIndex();
%             obj.displacementErrors_(iTimeStep,1) = time;
%             obj.displacementErrors_(iTimeStep,2) = meanErrorPerc;
%             obj.displacementErrors_(iTimeStep,3) = maxErrorPerc;
%             obj.displacementErrors_(iTimeStep,4) = minErrorPerc;
%         end
        
        function [obj, meanError, maxError, minError] = ...
                displacementAbsError...
                (...
                    obj,...
                    calculatedU,...
                    analyticU...
                )
            
            mesh = obj.mesh();
            N = mesh.numberOfElements;
            
%             residual = 0;
%             maxError = -inf;
%             minError = inf;
%             for iElem = 1:N
%                 nU = phi.internalField(iElem);    % numerical U
%                 aU = U.internalField(iElem);      % analytical U
%                 
%                 error = norm(nU - aU);
%                 
%                 if error > maxError
%                     maxError = error;
%                 end
%                 
%                 if error < minError
%                     minError = error;
%                 end
%                     
%                 residual = residual + error;
%             end
            
            dU = calculatedU.internalField() - analyticU.internalField();
            dUnorm = norm(dU);
            residual = sum(dUnorm.data);
            maxError = max(dUnorm.data);
            minError = min(dUnorm.data);
            euclideanNorm = fieldNorm(dU,2,true);
            
            disp('The absolute error = | U_calculated - U_analytic |:');
            
            meanError = residual/N;
            sError = sprintf('%0.2g', meanError);
            fprintf('Mean error U = %s\n', sError);
            
            sMaxError = sprintf('%0.2g', maxError);
            fprintf('Max  error U = %s\n', sMaxError);
            
            sMinError = sprintf('%0.2g', minError);
            fprintf('Min  error U = %s\n', sMinError);
            
            sEuclideanError = sprintf('%0.2g', euclideanNorm);
            fprintf('Euclidean  error U = %s\n', sEuclideanError);
            
            if euclideanNorm > 1e+3
                error('Simulation diverged: Euclidiean error is too high');
            end
        end
        
        function [obj, meanError, maxError, minError] = ...
                displacementNormalisedAbsError...
                (...
                    obj,...
                    calculatedU,...
                    analyticU...
                )
            %| U_calculated - U_analytic | / max(| U_analytic_field |)
            
            mesh = obj.mesh();
            N = mesh.numberOfElements;
            
            normField = analyticU.internalField().norm();
            maxDisp = fieldNorm(normField,Inf);
            
            if maxDisp == 0
                maxDisp = 1;
            end
            
%             if RUN_OPTIMISATIONS
%                 iPhi = phi.data; % Internal Phi data
%                 iU = U.data;     % Internal U data
%             end
%             
%             for iElem = 1:N
%                 
%                 if RUN_OPTIMISATIONS
%                     nU = iPhi(:,iElem);    % numerical U
%                     aU = iU(:,iElem);      % analytical U
%                 else
%                     nU = phi.internalField(iElem);    % numerical U
%                     aU = U.internalField(iElem);      % analytical U
%                 end
%                 
%                 error = norm(nU - aU)/maxDisp;
%                 
%                 if error > maxError
%                     maxError = error;
%                 end
%                 
%                 if error < minError
%                     minError = error;
%                 end
%                     
%                 residual = residual + error;
%             end
            
            error = norm(calculatedU.internalField() - analyticU.internalField())/maxDisp;
            residual = sum(error.data);
            maxError = max(error.data);
            minError = min(error.data);
            
            disp(['The normalised error = | U_calculated - U_analytic |', ...
                ' / max(| U_analytic_field |):']);
            
            meanError = residual/N;
            sError = sprintf('%0.2g', meanError);
            fprintf('Mean error U = %s\n', sError);
            
            sMaxError = sprintf('%0.2g', maxError);
            fprintf('Max  error U = %s\n', sMaxError);
            
            sMinError = sprintf('%0.2g', minError);
            fprintf('Min  error U = %s\n', sMinError);
        end
        
        function obj = displacementError_x(obj, phi, time)
            
            mesh = obj.mesh();
            x = obj.analyticalVol_x(time);
            N = mesh.numberOfElements;
            
            residual = 0;
            skippeds = 0;
            maxError = -inf;
            minError = inf;
            for iElem = 1:N
                nx = phi.internalField(iElem);    % numerical x
                ax = x.internalField(iElem);      % analytical x
                P = norm(ax);
                
                if P < 1e-5
                    skippeds = skippeds + 1;
                    continue;
                end
                
                error = norm(nx - ax) / P;
                
                if error > maxError
                    maxError = error;
                end
                
                if error < minError
                    minError = error;
                end
                    
                residual = residual + error;
            end
            sError = sprintf('%0.2g%%', residual/(N-skippeds)*100);
            fprintf('Mean error x = %s\n', sError);
            
            sMaxError = sprintf('%0.2g%%', maxError*100);
            fprintf('Max  error x = %s\n', sMaxError);
            
            sMinError = sprintf('%0.2g%%', minError*100);
            fprintf('Min  error x = %s\n', sMinError);
        end
        
        function obj = statistics(obj)
            
            time = obj.runTime().value();
            calculatedU = obj.solidModel().volU();
            analyticU = obj.analyticalVolU(time);
            
            [~, meanError, maxError, minError] = ...
                obj.displacementAbsError(calculatedU, analyticU);
            obj.displacementNormalisedAbsError(calculatedU, analyticU);
            
            mesh = obj.mesh();
            iTimeStep = mesh.timeMesh.timeIndex();
            if iTimeStep > 0
                obj.displacementErrors_(iTimeStep,1) = time;
                obj.displacementErrors_(iTimeStep,2) = meanError;
                obj.displacementErrors_(iTimeStep,3) = maxError;
                obj.displacementErrors_(iTimeStep,4) = minError;
            end
        end
        
        function obj = statistics_x(obj, solidModel)
            
            vol_x = solidModel.vol_x();
            obj.displacementError_x(vol_x, obj.runTime().value());
        end
        
        function printDisplacementErros(obj)
            
            for i=1:size(obj.displacementErrors_,1)
                
                time = obj.displacementErrors_(i,1);
                mean = obj.displacementErrors_(i,2);
                maxi = obj.displacementErrors_(i,3);
                mini = obj.displacementErrors_(i,4);
                
                fprintf(...
                    '%s  %s  %s  %s \n',...
                    sprintf('%0.2g', time),...
                    sprintf('%0.2g', mean),...
                    sprintf('%0.2g', maxi),...
                    sprintf('%0.2g', mini)...
                    )
            end
        end
        
        function solidModelHasEvolved(obj)
            %solidModelHasEvolved should be called after the solidModel has
            % evolved. Override this method as you wish.
            
            obj.statistics();
            
            obj.solidModelHasEvolved@Problem();
        end
    end
end

