% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.BlockCoupled;
%  params.Lambda = 1.45;
%  pathToMesh = '../meshes/unitCube3x3x3';
%  problem = kamojjala.Uniaxial(pathToMesh, params);
%  problem.solve();

checkRunningLive();

params = [];
params.classPrefix = 'kamojjala.Uniaxial';
params.runningOnConsole = false; 
params.numberOfTimeSteps = 10;
params.numberOfCorrectors = 2;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Cork216;
params.Lambda = 1.4;
params.useTraction = true;
%params.alternativeTolerance = 1e-3;
params.logFilename = [mfilename, '.log' ];

% 1 - {'unitCube16x16x16'}
% 2 - {'unitCube32x32x32'}
% 3 - {'unitCube3x3x3'}
% 4 - {'unitCube8x8x8'}
params.meshIdx = 3;

problem = caseRunner(params);

%notifySimulationEnded();
