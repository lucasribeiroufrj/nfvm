function volU = solveBlockCoupledLinearSystem(volU, A, B)

nTotalElements = size(A, 1)/3;
U = A\B;
volU = ...
    volU.setFromRawData ...
    ( ...
        reshape...
        (...
            U,...
            [3, nTotalElements] ...
        ),...
        1,...
        nTotalElements ...
    );

end