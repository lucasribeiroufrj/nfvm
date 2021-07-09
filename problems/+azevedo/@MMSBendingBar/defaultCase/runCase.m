% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModel = 'Segregated';
%  params.peakAmplitude = pi/7;
%  pathToMesh = '../meshes/unitCube3x3x3';
%  main(azevedo.MMSBendingBar(pathToMesh, params), params);

checkRunningLive();

params = [];
params.classPrefix = 'azevedo.MMSBendingBar';
params.runningOnConsole = false; 
params.peakAmplitude = pi*0.999;
params.useTraction = false;
params.numberOfTimeSteps = 3;
params.numberOfCorrectors = 1000;
params.solutionTolerance = 1e-6;
params.solidModelSolver = solidModelSolvers.List.Segregated;
params.materialLaw = materialLaws.List.NeoHookean;
params.material = materials.List.Steel;

if params.materialLaw == materialLaws.List.HookeanElastic ...
    || params.materialLaw == materialLaws.List.StVenant

    % Small deformation/strain.
    params.peakAmplitude = pi/9;
end

% 1 - {'unitCube3x3x3'}
% 2 - {'unitCube8x8x8'}
% 3 - {'unitCubeAnsys300'}
params.meshIdx = 2;

% Callback function called after problem object has received the
% message informing about solidModel evolution/update. We use to stop the
% simulation when simulation time is endTime/2 since the simulation is
% periodic.
params.solidModelhasEvolvedCallback = @solidModelhasEvolvedCallback;

params.logFilename = [mfilename, '.log' ];

problem = caseRunner(params);

problem.printDisplacementErros();

%notifySimulationEnded();

function solidModelhasEvolvedCallback(theProblem)

    if theProblem.runTime.value() >= theProblem.runTime.endTime()/2
        theProblem.runTime.setEndTime(0);
    end
end
