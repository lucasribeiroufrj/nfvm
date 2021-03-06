function mechanicalModel = ...
    createMechanicalModel(materialLaw, solidModelSolver, varargin)
%createMechanicalModel Summary of this function goes here
%   Detailed explanation goes here

if solidModelSolver == solidModelSolvers.List.NonLinearBlockCoupled
    
    mechanicalModel = ...
        mechanicalModels.NonLinearBlockCoupledMechanicalModel...
        (...
            materialLaw,...
            varargin{:}...
        );
    
elseif solidModelSolver == solidModelSolvers.List.BlockCoupled
    
    mechanicalModel = ...
        mechanicalModels.BlockCoupledMechanicalModel...
        (...
            materialLaw,...
            varargin{:}...
        );
else
    
    mechanicalModel = ...
        mechanicalModels.SegregatedMechanicalModel...
        (...
            materialLaw,...
            varargin{:}...
        );
end

end

