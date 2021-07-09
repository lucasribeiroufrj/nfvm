function [ newPhi ] = linearSolverDirect( ~, A, B, ~ )
%linearSolverDirect Summary of this function goes here
%   Direct: A * phi = B  (x A^-1 both sides)

newPhi = (A\B).';

end

