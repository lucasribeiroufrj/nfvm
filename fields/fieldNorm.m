function residual = fieldNorm(genericField, p, ignoreBoundary)
%fieldNorm of a geometric field.
%   fieldNorm(genericField, p) returns the Frobenius norm of the 
% genericField if p = 'fro'. Other values for p might be used
% as in MATLAB's norm function.
%   genericField: (vol/surface)TensorField, (vol/surface)VectorField.
%

if nargin == 1
    p = 2;
    ignoreBoundary = false;
elseif nargin == 2
    ignoreBoundary = false;
end

if isa(genericField, 'volTensorField') ...
        || isa(genericField, 'surfaceTensorField')

    residual = tensorFieldNorm(genericField, p);

    if ~ignoreBoundary
        
        for iBoundary = 1:genericField.mesh().numberOfBoundaries()

            boundary = genericField.boundaryField(iBoundary);
            residual = residual + tensorFieldNorm(boundary, p);
        end
    end
elseif isa(genericField, 'volVectorField') ...
        || isa(genericField, 'surfaceVectorField')

    residual = vectorFieldNorm(genericField, p);

    if ~ignoreBoundary
        
        for iBoundary = 1:genericField.mesh().numberOfBoundaries()
            
            boundary = genericField.boundaryField(iBoundary);
            residual = residual + vectorFieldNorm(boundary, p);
        end
    end
    
elseif isa(genericField, 'tensorField') ...
    
    residual = tensorFieldNorm(genericField, p);
    
elseif isa(genericField, 'vectorField') ...
    
    residual = vectorFieldNorm(genericField, p);
    
elseif isa(genericField, 'scalarField') ...
    
    residual = scalarFieldNorm(genericField, p);
else
    error('Not implemented');
end

end


