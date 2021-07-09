% Use the script this script or type the command below to run a case:
%
% For example:
%  params = [];
%  params.solidModel = 'Segregated';
%  params.Lambda = 1.45;
%  pathToMesh = '../../meshes/unitCube3x3';
%  main(kamojjala.Uniaxial2d(pathToMesh, params), params);
%%

checkRunningLive();

solidModels = {'Segregated', 'NonLinearBlockCoupled'};
laws = [materialLaws.List.NeoHookean];
meshes = [6, 11, 1, 4, 8];

for s = 1:length(solidModels)
    for l = 1:length(laws)
        for i=1:length(meshes)

            try
            params = [];
            params.classPrefix = 'kamojjala.Uniaxial2d';
            params.runningOnConsole = false; 

            params.solidModel = solidModels{s};

            params.Lambda = 2.0;
            params.useTraction = false;
            params.numberOfTimeSteps = 1;
            params.numberOfCorrectors = 20000;
            params.useSymmetryBoundaries = false;

            params.law = laws(l);

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
            params.meshIdx = meshes(i);
            
            fprintf('Simulation params:\n');
            params %#ok<NOPTS>

            caseRunner(params);

            notifySimulationEnded();

            catch msg

                errMsg = ['Error => ', msg.message];
                notifySimulationEnded({'message', errMsg})
                disp(errMsg);
            end
        end
    end
end
%%
