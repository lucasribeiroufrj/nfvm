% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  pathToMesh = '../meshes/plateHole250cells'
%  problem = demirdzic.PlateHole(pathToMesh, params);
%  problem.solve();

params = [];
params.classPrefix = 'demirdzic.PlateHole';
params.runningOnConsole = false;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.HookeanElastic;
params.material = materials.List.Rubber;
params.logFilename = [mfilename, '.log' ];
params.meshIdx = 1;

% 1 - 'plateHole1000cells'
% 2 - 'plateHole10cells'
% 3 - 'plateHole2250cells'
% 4 - 'plateHole250cells'
% 5 - 'plateHole4000cells'
% 6 - 'plateHoleTwoblocks10x10'
% 7 - 'plateHoleTwoblocks20x20'

problem = caseRunner(params);
