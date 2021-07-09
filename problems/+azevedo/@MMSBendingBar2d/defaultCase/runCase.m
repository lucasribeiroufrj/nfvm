% Use this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  params.peakAmplitude = pi/9;
%  pathToMesh = '../meshes/unitCube3x3';
%  main(azevedo.MMSBendingBar2d(pathToMesh, params), params)
tic
checkRunningLive();

params = [];
params.classPrefix = 'azevedo.MMSBendingBar2d';
params.runningOnConsole = false;
params.solidModelSolver = solidModelSolvers.List.NonLinearBlockCoupled;
params.peakAmplitude = pi/7;
params.numberOfTimeSteps = 3;%20;
params.numberOfCorrectors = 1000;
params.useTraction = false;
params.alternativeTolerance = 1e-7;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Cork216;
%params.alternativeTolerance = 1e-3;
params.logFilename = [mfilename, '.log' ];

% 1 - {'bendingBar8x64'}
% 2 - {'bendingBar8x16'}
% 3 - {'unitCube16x16'}
% 4 - {'unitCube32x32'}
% 5 - {'unitCube3x3'}
% 6 - {'unitCube64x64'}
% 7 - {'unitCube80x80'} 
% 8 - {'unitCube8x8'}
params.meshIdx = 8;

if params.meshIdx == 1
    
    params.barHeight = 8;
    
elseif params.meshIdx == 2
    
    params.barHeight = 2;
end

problem = caseRunner(params);
toc
%notifySimulationEnded();