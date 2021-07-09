% Use this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  pathToMesh = '../meshes/cantileverBeam60x3x1';
%  main(cardiff.CantileverBeam2d(pathToMesh, params), params);

checkRunningLive();

params = [];
params.classPrefix = 'cardiff.CantileverBeam2d';
params.runningOnConsole = false;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Steel;
params.numberOfCorrectors = 100000;
%params.alternativeTolerance = 1e-7;
params.logFilename = [mfilename, '.log' ];

if params.solidModelSolver == solidModelSolvers.List.NonLinearBlockCoupled
    params.numberOfCorrectors = 1;
end

% 1 - {'cantileverBeam100x5x1'}
% 2 - {'cantileverBeam300x15x1'}
% 3 - {'cantileverBeam3x3x1'}
% 4 - {'cantileverBeam60x3x1'}
params.meshIdx = 1;

problem = caseRunner(params);

%notifySimulationEnded();