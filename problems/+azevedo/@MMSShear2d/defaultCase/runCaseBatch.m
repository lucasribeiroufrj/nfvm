% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModel = 'Segregated';
%  params.Lambda = 1.3;
%  pathToMesh = '../meshes/unitCube3x3';
%  main(azevedo.MMSShear2d(pathToMesh, params), params);

checkRunningLive();

meshNames = {'unitCube3x3','unitCube8x8','unitCube16x16',...
    'unitCube32x32','unitCube64x64'};
useTraction = {false, true};
solidModels = {'Segregated', 'NonLinearBlockCoupled'};
laws = [MaterialLaws.NeoHookean];
meshes = [3, 5, 1, 2, 4];

for t = 1:length(useTraction)
    for s = 1:length(solidModels)
        for l = 1:length(laws)
            for i=1:length(meshes)
                try
                    params = [];
                    params.classPrefix = 'azevedo.MMSShear2d';
                    params.runningOnConsole = false; 

                    params.solidModel = solidModels{s};

                    params.shearFactor = 0.45;
                    params.useTraction = useTraction{t};
                    params.numberOfTimeSteps = 1;
                    params.numberOfCorrectors = 20000;

                    params.law = laws(l);

                    % 1 - {'unitCube16x16'}
                    % 2 - {'unitCube32x32'}
                    % 3 - {'unitCube3x3'}
                    % 4 - {'unitCube64x64'}
                    % 5 - {'unitCube8x8'}
                    params.meshIdx = meshes(i);

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

                figName = sprintf('UseTraction_%d_%s%s_Mesh_%s', ...
                    params.useTraction,...
                    params.solidModel,...
                    params.law,...
                    meshNames{i});
                savefig(figName);
            end
        end
    end
end
