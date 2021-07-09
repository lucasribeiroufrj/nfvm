function value = scalarFieldNorm(field, p)
%scalarFieldNorm norm of a vectorField.
%   scalarFieldNorm(field, p) returns the p norm of the 
% just as MATLAB's norm function.
% Ex.: 
%   if p =  Inf => residual = max(abs(field))
%   if p = -Inf => residual = min(abs(field))
%   field: a vectorField object.

value = norm(field.data, p);

end
