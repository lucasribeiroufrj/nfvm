function checkDiagonalDominance( A )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if all((2*abs(diag(A))) >= sum(abs(A),2))
    fprintf('The matrix is diagonal dominant.\n');
else
    fprintf('The matrix is not diagonal dominant\n');
end

end

