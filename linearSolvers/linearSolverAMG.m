function [ newPhi ] = linearSolverAMG(phi, A, B, config)
%linearSolverAMG Summary of this function goes here
%   Detailed explanation goes here

newPhi = phi*0; % Allocate space

opt = amgset;                               % Default option file.
opt = amgset(opt,'coarsest',2);             % Det the number of levels.
opt = amgset(opt,'PreCond','pcg');          % Det the Krylov method.
opt = amgset(opt,'PrintOnScreen','off');    % Turn off the option to print
                                            % the log on the screen.
opt = amgset(opt,'SaveCsn','off');          % Save the set of coarse-grid 
                                            % points.
opt = amgset(opt,'CsnType','amg');          % Choose the coarsening method.
opt = amgset(opt,'TolAMG', config.linearSystemTolerance);
opt = amgset(opt, 'RelPara', config.amgRelaxation);

for cmp=1:3
    newPhi(cmp,:) = ...
        amg(A,phi(cmp,:).',B(:,cmp),opt);
end

end

