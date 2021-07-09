function [ newPhi ] = linearSolverICCG(phi, A, B, config)
%linearSolverICCG Summary of this function goes here
%   Detailed explanation goes here
%
%  TODO:
% - We are converting the input matrix to a sparse matrix. This is very
% slow!
%

newPhi = phi*0; % Allocate space
sparseA = sparse(-A);
converged = true;

for cmp=1:3

    L = ichol(sparseA,struct('michol','on'));
    [x,flag,~,iter,relResVec] = ...
        pcg...
        (...
            sparseA, -B(:,cmp), config.linearSystemTolerance, 100, L, L', ...
            phi(cmp,:).'...
        );

    newPhi(cmp,:) = x;

    if flag == 0
        converged = converged & true;
    else
        converged = converged & false;
    end

    if false
    if ~RUNNING_ON_CONSOLE && false
        figure(9);
        clf(9)
        semilogy(0:iter,relResVec(1:iter+1)/norm(B(:,cmp)),'k.');
        legend(sprintf('Component x_%d', cmp));
        xlabel('iteration number');
        ylabel('relative residual');
    end
    end
end

msg = 'ICCG: Inner iteration ';
if converged
    disp([ msg 'converged.' ]);
else
    disp([ msg 'did not converge.' ]);
end

end

