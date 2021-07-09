function config = configGenerator()
%CONFIGGENERATOR Summary of this function goes here
%   Detailed explanation goes here

% Heavly dependent on the mesh. Two iterations seems to avoid 
% oscillations [Pag. 277 - Moukalled].
config.maxGradIterComputation = 2;

% Maximum number of momentum correctors
config.nCorrectors=150;

% For incremental approaches
config.nDisplacementIncrements=1;

%config.linearSolver = 'Direct';
config.linearSolver = 'ICCG';
%config.linearSolver = 'AMG';

%config.gradInterpolationScheme = 'Least-Square';
config.gradInterpolationScheme = 'Gauss-Green';

% Linear system relaxation
config.relaxation = 1.0;

% Solution tolerance for displacement
config.solutionTolerance = 1e-6;
config.alternativeTolerance = 1e-7;

% Solution tolerance linear system
config.linearSystemTolerance = 1e-8;

config.printFrequency = 1;
config.applyUnderRelaxation = false;
config.amgRelaxation = 1.0;

end

