% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  params.deltaLength = 1.3;
%  pathToMesh = '../meshes/unitCube3x3/';
%  main(azevedo.Uniaxial2d(pathToMesh, params), params);

checkRunningLive();

params = [];
params.classPrefix = 'azevedo.Uniaxial2d';
params.runningOnConsole = false; 
params.numberOfTimeSteps = 1;
params.numberOfCorrectors = 50;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Cork216;
params.deltaLength = 1.0;
%params.alternativeTolerance = 1e-3;
params.logFilename = [mfilename, '.log' ];

% 1 - {'unitCube3x3'}
% 2 - {'unitCube8x8'}
params.meshIdx = 2;

if (params.solidModelSolver == solidModelSolvers.List.NonLinearBlockCoupled)
    params.numberOfCorrectors = 1;
end

problem = caseRunner(params);

%notifySimulationEnded();
