function residual = tensorFieldNorm(field, p)
%tensorFieldNorm norm of a tensorField.
%   tensorFieldNorm(field, p) returns the Frobenius norm of the 
% field if p = 'fro'. Other values for p might be used
% as in MATLAB's norm function.
%   field: a tensorField object.

data = field.data;
residual = 0;

for iElement = 1:field.numberOfElements()
    residual = residual + norm(data(:,:,iElement), p);
end

end
