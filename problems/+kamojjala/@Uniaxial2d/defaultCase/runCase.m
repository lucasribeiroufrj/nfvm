% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.BlockCoupled;
%  params.Lambda = 1.45;
%  pathToMesh = '../meshes/unitCube3x3';
%  main(kamojjala.Uniaxial2d(pathToMesh, params), params);
%%

checkRunningLive();

try
    params = [];
    params.classPrefix = 'kamojjala.Uniaxial2d';
    params.runningOnConsole = false; 

    params.solidModelSolver = solidModelSolvers.List.NonLinearBlockCoupled;

    params.Lambda = 0.85;
    params.useTraction = true;
    params.numberOfTimeSteps = 1;
    params.numberOfCorrectors = 1000;
    params.useSymmetryBoundaries = false; % symmetry causes "skids".

    params.materialLaw = materialLaws.List.NeoHookean;
    params.material = materials.List.Cork216;

    % 1 - {'unitCube16x16'}
    % 2 - {'unitCube16x16_grading'} 
    % 3 - {'unitCube17x17'}
    % 4 - {'unitCube32x32'}
    % 5 - {'unitCube33x33'}
    % 6 - {'unitCube3x3'}
    % 7 - {'unitCube5x5'}
    % 8 - {'unitCube64x64'}
    % 9 - {'unitCube65x65'}
    % 10 - {'unitCube80x80'}
    % 11 - {'unitCube8x8'}
    % 12 - {'unitCube8x8_grading'}
    % 13 - {'unitCube9x9'}
    params.meshIdx = 1;
    
    params.logFilename = [mfilename, '.log' ];

    fprintf('Simulation params:\n');
    params %#ok<NOPTS>
    
    if (params.solidModelSolver == solidModelSolvers.List.NonLinearBlockCoupled)
        params.numberOfCorrectors = 10;
    end

    problem = caseRunner(params);

    notifySimulationEnded();

catch msg

    errMsg = ['Error => ', msg.message];
    disp(msg.stack(1));
    notifySimulationEnded({'message', errMsg})
    disp(errMsg);
end
