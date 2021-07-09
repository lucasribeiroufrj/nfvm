% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  params.Lambda = 1.3;
%  pathToMesh = '../meshes/unitCube3x3';
%  main(azevedo.MMSSpring2d(pathToMesh, params), params);

checkRunningLive();

params = [];
params.classPrefix = 'azevedo.MMSSpring2d';
params.runningOnConsole = false; 
params.numberOfTimeSteps = 1;
params.numberOfCorrectors = 1000;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Cork216;
%params.alternativeTolerance = 1e-3;
params.logFilename = [mfilename, '.log' ];

% 1 - {'unitCube16x16'}
% 2 - {'unitCube32x32'}
% 3 - {'unitCube3x3'}
% 4 - {'unitCube8x8'}
params.meshIdx = 1;

if params.solidModelSolver == solidModelSolvers.List.NonLinearBlockCoupled
    
    params.numberOfCorrectors = 1;
end

problem = caseRunner(params);

%notifySimulationEnded();
