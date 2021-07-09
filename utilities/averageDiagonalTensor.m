function average = averageDiagonalTensor(A)
% averageDiagonalTensor(A) returns the average of the diagonal tensor
%   elements of the block coupled matrix A.
    
    dofs = size(A,1);
    
    if mod(dofs,3) ~= 0
        
        ME = MException('averageDiagonalTensor:invalidInput',...
        'The input matrix must be multiple of 3');
        throw(ME)
    end
    
    nTotalElements = dofs/3;
    average = zeros(3,3);
    
    for iElement = 1:nTotalElements
        
        i = getTensorialIndex(iElement);
        average = average + A(i,i);
    end
    
    average = average/nTotalElements;
end
