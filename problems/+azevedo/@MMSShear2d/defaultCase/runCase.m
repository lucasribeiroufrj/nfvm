% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  params.Lambda = 1.3;
%  pathToMesh = '../meshes/unitCube3x3';
%  problem = azevedo.MMSShear2d(pathToMesh, params);
%  problem.solve(params);

checkRunningLive();

params = [];
params.classPrefix = 'azevedo.MMSShear2d';
params.runningOnConsole = false; 
params.numberOfTimeSteps = 1;
params.numberOfCorrectors = 1000;
params.solidModelSolver = solidModelSolvers.List.NonLinearBlockCoupled;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Cork216;
params.shearFactor = 0.45;
params.useTraction = true;
%params.alternativeTolerance = 1e-3;
params.logFilename = [mfilename, '.log' ];

% 1 - {'unitCube16x16'}
% 2 - {'unitCube32x32'}
% 3 - {'unitCube3x3'}
% 4 - {'unitCube64x64'}
% 5 - {'unitCube8x8'}
params.meshIdx = 4;

problem = caseRunner(params);

notifySimulationEnded();
