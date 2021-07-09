% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  params.Lambda = 1.3;
%  pathToMesh = '../meshes/unitCube3x3x3';
%  main(azevedo.Uniaxial(pathToMesh, params), params);

checkRunningLive();

params = [];
params.classPrefix = 'azevedo.Uniaxial';
params.runningOnConsole = false; 
params.deltaLength = 0.85;
params.numberOfTimeSteps = 1;
params.numberOfCorrectors = 50;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Cork216;

% 1 - {'unitCube3x3x3'}
% 2 - {'unitCube8x8x8'}
% 3 - {'unitCubeAnsys300'}
params.meshIdx = 2;

params.logFilename = [mfilename, '.log' ];

if (params.solidModelSolver == solidModelSolvers.List.NonLinearBlockCoupled)
    params.numberOfCorrectors = 1;
end

problem = caseRunner(params);

%notifySimulationEnded();
