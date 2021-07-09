function newPhi = solveLinearSystem(phi, A, B, config)
%solveEquation Solve the linear system
%   Detailed explanation goes here

    % Patankar's implicit method [p. 150, Moukalled - 2015; 
    % eq. 3.56, Muzaferija - 1994]
    function applyUnderRelaxation()
        
        alpha = config.relaxation;
        
        for iElem=1:size(A,1)
            
            Ac = A(iElem, iElem) / alpha;
            A(iElem, iElem) = Ac;
            B(iElem,:) = B(iElem,:) + (1-alpha)*Ac*phi(:,iElem).'; 
        end
    end

if config.applyUnderRelaxation
    applyUnderRelaxation();
end

if strcmp(config.linearSolver,'ICCG')
    
    newPhi = linearSolverICCG(phi, A, B, config);
    
elseif strcmp(config.linearSolver, 'AMG')
    
    newPhi = linearSolverAMG(phi, A, B, config);
    
elseif strcmp(config.linearSolver, 'Direct')
    
    newPhi = linearSolverDirect(phi, A, B, config);
else
    ME = MException('Error: The linear solver %s is not supported.', ...
        config.linearSolver);
    throw(ME);
end

if ~isreal(newPhi(1,1))
    error('Imaginary solution detected.');
end

end
