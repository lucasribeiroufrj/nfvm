% Use this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModelSolver = solidModelSolvers.List.Segregated;
%  params.peakAmplitude = pi/9;
%  pathToMesh = '../meshes/unitCube3x3';
%  main(azevedo.MMSBendingBar2d(pathToMesh, params), params)

checkRunningLive();

meshNames = {'unitCube3x3','unitCube8x8','unitCube16x16',...
    'unitCube32x32','unitCube64x64'};
useTraction = {false, true};
solidModelSolvers = ...
    {...
        solidModelSolvers.List.Segregated...
        solidModelSolvers.List.NonLinearBlockCoupled...
    };
materialLaws = [materialLaws.List.NeoHookean];
meshes = [5, 8, 3, 4, 6];
params.material = materials.List.Cork216;

for t = 1:length(useTraction)
    for s = 1:length(solidModelSolvers)
        for l = 1:length(materialLaws)
            for i=1:length(meshes)
                try
                    params = [];
                    params.classPrefix = 'azevedo.MMSBendingBar2d';
                    params.runningOnConsole = false; 

                    params.solidModelSolver = solidModelSolvers{s};

                    params.peakAmplitude = pi/7;
                    params.numberOfTimeSteps = 3;
                    params.numberOfCorrectors = 10000;
                    params.useTraction = useTraction{t};
                    params.alternativeTolerance = 1e-7;

                    params.materialLaw = materialLaws(l);

                    % 1 - {'bendingBar8x64'}
                    % 2 - {'cubeBending8x8'}
                    % 3 - {'unitCube16x16'}
                    % 4 - {'unitCube32x32'}
                    % 5 - {'unitCube3x3'}
                    % 6 - {'unitCube64x64'}
                    % 7 - {'unitCube80x80'} 
                    % 8 - {'unitCube8x8'}
                    params.meshIdx = meshes(i);

                    if params.meshIdx == 1

                        params.barHeight = 8;

                    elseif params.meshIdx == 2

                        params.barHeight = 2;
                    end
                    
                    % To print information about each run
                    params.beforeCallMainCallback = @() disp(params);

                    params.logFilename = [mfilename, '.log' ];

                    problem = caseRunner(params);

                    notifySimulationEnded();
                    
                catch msg

                    errMsg = ['Error => ', msg.message];
                    notifySimulationEnded({'message', errMsg})
                    disp(errMsg);
                end

                figName = sprintf('traction_%d_%s%s_Mesh_%s_ts_%d', ...
                    params.useTraction,...
                    params.solidModelSolver,...
                    params.materialLaw,...
                    meshNames{i},...
                    params.numberOfTimeSteps);
                savefig(figName);
            end
        end
    end
end
