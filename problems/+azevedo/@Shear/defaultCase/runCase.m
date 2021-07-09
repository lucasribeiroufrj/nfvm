% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  params.deltaLength = 0.45;
%  pathToMesh = '../meshes/unitCube3x3x3';
%  main(azevedo.Shear(pathToMesh, params), params);

checkRunningLive();

params = [];
params.classPrefix = 'azevedo.Shear';
params.runningOnConsole = false; 
params.numberOfTimeSteps = 1;
params.numberOfCorrectors = 1000;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Cork216;
params.deltaLength = 0.45;
%params.alternativeTolerance = 1e-3;
params.logFilename = [mfilename, '.log' ];

% 1 - {'unitCube3x3x3'}
% 2 - {'unitCube8x8x8'}
params.meshIdx = 2;

if params.solidModelSolver == solidModelSolvers.List.NonLinearBlockCoupled
    
    params.numberOfCorrectors = 1;
end

problem = caseRunner(params);

notifySimulationEnded();