function residual = vectorFieldNorm(field, p)
%vectorFieldNorm norm of a vectorField.
%   vectorFieldNorm(field, p) returns the Frobenius norm of the 
% field if p = 'fro'. Other values for p might be used
% as in MATLAB's norm function.
%   field: a vectorField object.

data = field.data;
residual = 0;

for iElement = 1:field.numberOfElements()
    residual = residual + norm(data(:,iElement), p);
end

end

