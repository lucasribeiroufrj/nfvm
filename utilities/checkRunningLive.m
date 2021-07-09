function checkRunningLive()
%CHECKRUNNINGLIVE Call it to avoid running live
%   Detailed explanation goes here

stackTrace = dbstack(1);%, '-completenames');

if strcmp(stackTrace(1).name, 'LiveEditorEvaluationHelperESectionEval')
    error('Sorry, cannot run live. Use Run buttom.');
end

end

