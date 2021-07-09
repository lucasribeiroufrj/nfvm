% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  pathToMesh = '../meshes/unitCube3x3x3';
%  problem = demirdzic.FreeCube(pathToMesh, params);
%  problem.solve();

checkRunningLive();

params = [];
params.classPrefix = 'demirdzic.FreeCube';
params.runningOnConsole = false;
params.solidModelSolver = solidModelSolvers.List.NonLinearBlockCoupled;
params.materialLaw = materialLaws.List.DuhamelNeumann;
params.material = materials.List.Rubber;
params.numberOfCorrectors = 1000;
%params.alternativeTolerance = 1e-9;
params.meshIdx = 1;
params.logFilename = [mfilename, '.log' ];

problem = caseRunner(params);

%notifySimulationEnded();
