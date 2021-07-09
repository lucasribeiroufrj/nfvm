function notifySimulationEnded(params)
%notifySimulationEnded Send notification to OS.
%   Send msg to notification center of the OS.

[status, result] = system('uname');
name = result(1:end-1);

if status == 0 

    if nargin == 0
        params = [];
    end
    
    message = ...
        tryGetOrDefault(params, 'message', 'Simulation has ended.');
    
    title = ...
        tryGetOrDefault(params, 'title', 'MATLAB');
    
    if strcmp(name, 'Darwin')
        
        system(sprintf(['osascript -e "display notification '...
            '\\"%s\\" with ' ...
            'title \\"%s\\""'], message, title));
    else
        fprintf('Notification center for %s OS is not implemented\n.',...
            name);
    end 
else
    disp('The name of the operational system could not be retrieved');
end

end

