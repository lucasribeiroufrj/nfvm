function printFooter()
%PRINTFOOTER Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');
fprintf('Simulation ended\n');
fprintf('Date     : %s\n', date);
fprintf('Time     : %s\n', datestr(now, 'HH:MM:SS'));
fprintf('\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');

end

