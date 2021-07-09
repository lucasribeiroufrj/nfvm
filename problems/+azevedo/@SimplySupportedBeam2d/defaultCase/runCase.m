% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  pathToMesh = '../meshes/simplySupportedBeam60x3';
%  main(azevedo.SimplySupportedBeam2d(pathToMesh, params), params);

checkRunningLive();

params = [];
params.classPrefix = 'azevedo.SimplySupportedBeam2d';
params.runningOnConsole = false; 
params.numberOfTimeSteps = 1;%20;
params.numberOfCorrectors = 50000;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Rubber;
params.meshIdx = 2;
%params.alternativeTolerance = 1e-3;
params.logFilename = [mfilename, '.log' ];

if params.solidModelSolver == solidModelSolvers.List.NonLinearBlockCoupled
    
    params.numberOfCorrectors = 1;
end

problem = caseRunner(params);

%notifySimulationEnded();
