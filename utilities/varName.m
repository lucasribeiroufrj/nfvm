function [ out ] = varName( ~ )
%varName Get the name of the variable
%   Then try for example:
%
%   age = 18;
%   name = varname(age)
%   disp(name)

out = inputname(1);

end

