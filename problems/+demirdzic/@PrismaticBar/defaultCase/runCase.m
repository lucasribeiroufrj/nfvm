% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  pathToMesh = '../meshes/3x3x25';
%  problem = demirdzic.PrismaticBar(pathToMesh, params);
%  problem.solve();

checkRunningLive();

params = [];
params.classPrefix = 'demirdzic.PrismaticBar';
params.runningOnConsole = false; 
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.HookeanElastic;
params.material = materials.List.Steel;
params.useTraction = true;
%params.alternativeTolerance = 1e-3;
params.logFilename = [mfilename, '.log' ];

params.meshIdx = 1;

problem = caseRunner(params);

%notifySimulationEnded();
