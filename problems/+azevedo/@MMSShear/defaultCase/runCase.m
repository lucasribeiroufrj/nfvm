% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  params.shearFactor = 1.3;
%  pathToMesh = '../meshes/unitCube3x3x3';
%  problem = azevedo.MMSShear(pathToMesh, params);
%  problem.solve();

checkRunningLive();

params = [];
params.classPrefix = 'azevedo.MMSShear';
params.runningOnConsole = false;
params.numberOfTimeSteps = 1;
params.numberOfCorrectors = 1000;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.HookeanElastic;
params.material = materials.List.Cork216;
params.shearFactor = 0.45;
params.useTraction = true;
%params.alternativeTolerance = 1e-3;
params.logFilename = [mfilename, '.log' ];

% 1 - {'unitCube3x3x3'}
% 2 - {'unitCube8x8x8'}
% 3 - {'unitCubeAnsys300'}
params.meshIdx = 2;

problem = caseRunner(params);

%notifySimulationEnded();
