% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  pathToMesh = '../meshes/unitCube3x3x3';
%  main(demirdzic.ConstrainedCube(pathToMesh, params), params);

params = [];
params.classPrefix = 'demirdzic.ConstrainedCube';
params.runningOnConsole = false;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.DuhamelNeumann;
params.material = materials.List.Rubber;
%params.numberOfCorrectors = 2;
params.meshIdx = 1;

problem = caseRunner(params);
