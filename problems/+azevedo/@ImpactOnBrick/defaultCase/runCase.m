
% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModel = 'Segregated';
%  pathToMesh = '../meshes/brick4x4x4';
%  main(azevedo.ImpactOnBrick(pathToMesh, params), params);

params = [];
params.classPrefix = 'azevedo.ImpactOnBrick';
params.runningOnConsole = false;
params.numberOfCorrectors = 50000;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Cork216;
params.alternativeTolerance = 1e-7;

% 1 - {'brick4x4x4'}
% 2 - {'brick4x4x4Squared'}
% 3 - {'brick5x5x5'}
% 4 - {'brick8x8x8'}
% 5 - {'brickFromAnsys1340'}
% 6 - {'brickSym18x10x18'}
% 7 - {'brickSym4x4x4'}
% 8 - {'brickSym4x4x4Squared'}
% 9 - {'brickSym8x8x8'}
% 10- {'brickSymAnsys4682'}

params.meshIdx = 4;
params.logFilename = [mfilename, '.log' ];

problem = caseRunner(params);

%notifySimulationEnded();