function [hasConverged, residual] = ...
    computeResiduals...
    (...
        iCorrector,...
        volPhi,...
        volPhi_prevIter,...
        volPhi_o,...
        linearSolverInitialResidual,...
        config...
    )
%computeResiduals Summary of this function goes here
%   Detailed explanation goes here

solutionTolerance = config.solutionTolerance;
alternativeTolerance = config.alternativeTolerance;

    function value = vecnorm(M, p, direction)
        % vecnorm Vector-wise norm
        %   p = 2 -> Euclidean norm
        %   direction = 1 -> Column norm
        
        value = sqrt(sum(abs(M).^p, direction));
    end

    function value = euclideanNorm(vectorArray)
        
        value = vecnorm(vectorArray,2,1);
    end

    function [ hasConverged, residualDeltaD ] = cardiffConvergence()
        %cardiffConvergence This function comes from:
        %   solids4FoamModels/solidModels/solidModel/solidModel.C
        %
        % TODO
        % - this function could be in a separated file
        % - SMALL as a global value
        
        SMALL = 1.0e-15;
        
        phi = volPhi.internalField().data;
        phi_o = volPhi_o.internalField().data;
        phi_prevIter = volPhi_prevIter.internalField().data;
        
        denominator = max(euclideanNorm(phi - phi_o));
        
        if denominator < SMALL
            
            denominator = max(max(euclideanNorm(phi)), SMALL);
        end
        
        normArray = euclideanNorm(phi - phi_prevIter);
        
        residualDeltaD = max( normArray ) / denominator;
        
        if     linearSolverInitialResidual <  solutionTolerance...
            && residualDeltaD < solutionTolerance
            
            disp('=== Both residual have converged.')
            hasConverged = true;
            
        elseif residualDeltaD < alternativeTolerance
            
            disp('=== The relative residual has converged.')
            hasConverged = true;
            
        elseif linearSolverInitialResidual < alternativeTolerance
            
            disp('=== The linear solver residual has converged.')
            hasConverged = true;
        else
            hasConverged = false;
        end
        
        % TODO
        linearSolverIter = -1;

        oIter = sprintf('%d', iCorrector);
        oIterMax = sprintf('%d', config.nCorrectors);
        lsitRes = sprintf('%d', linearSolverInitialResidual);
        rDD = sprintf('%d', residualDeltaD);
        lsIter = sprintf('%d', linearSolverIter);
        disp('=== OuterIter, residual, relativeResidual, linearSolverIter')
        fprintf('=== %s/%s, %s, %s, %s\n', ...
            oIter, oIterMax, lsitRes, rDD, lsIter);

        solTor = sprintf('%0.0e', solutionTolerance);
        fprintf('=== Solution Tolerance: %s\n', solTor);

        altTor = sprintf('%0.0e', alternativeTolerance);
        fprintf('=== Alternative Tolerance: %s\n', altTor);
    end

fprintf('========================================================\n');

[ hasConverged, residual ] = cardiffConvergence();

fprintf('========================================================\n');

end

