function problem = caseRunner(classPrefix, params, meshIdx)
% caseRunner: helper function used to run built-in cases.
%   caseRunner(classPrefix, params, meshIdx):
%       ex.:
%       classPrefix = 'azevedo.ImpactOnBrick';
%       params = struct('runningOnConsole', 0);
%       meshIdx = 1;
%
%   Obs.: Prefer caseRunner(params) 
%   and for example, 
%       params.classPrefix = 'azevedo.ImpactOnBrick'
%       params.meshIdx = 1
%
%   For more information check the script run inside a built-in problem
%   folder.

%% BEG -  Get the available  meshes, then choose one.
if nargin == 1
    params = classPrefix;
    meshIdx = tryGetOrDefault(params, 'meshIdx', 1);
    classPrefix = params.classPrefix;
end
meshes = Problem.getNameOfCases(classPrefix);
mesh = meshes{meshIdx};
%% END

% Mount the mesh path.
root = Problem.getRootFolder(classPrefix);
fullPath = [root filesep 'meshes' filesep mesh];

% Create a problem.
problem = eval([ classPrefix '(''' fullPath ''', params)']);

% Log run
diary(tryGetOrDefault(params, 'logFilename', 'log.nFVM'));

% Set up a callback
beforeCallMainCallback = ...
    tryGetOrDefault(params, 'beforeCallMainCallback', @()[]);
beforeCallMainCallback();

% Run
problem.solve();

% Set up a callback
afterCallMainCallback = ...
    tryGetOrDefault(params, 'afterCallMainCallback', @()[]);
afterCallMainCallback();

% Touggle off the logging
diary

end
