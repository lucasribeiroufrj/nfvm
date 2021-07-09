function [LHS, RHS] = fvm_d2dt2(density, phi_o, phi_oo, problem)
%fvm_d2dt2 For each finite volume computes the temporal term contribution
%   to the final linear system. This is a first-order Euler implicit  using
%   the current and two previous time-step values. Using talyor:
%       phi_o = phi - deltaT*ddt(phi) + (deltaT)^2/2*d2dt2(phi)
%                + O(deltaT^3);                                     (I)
%       phi_oo = phi - (deltaT+deltaT_o)*ddt(phi) 
%                + (deltaT+deltaT_o)^2/2*d2dt2(phi) + O(deltaT^3).  (II)
%
%   Such that: 
%       phi_o  = phi(t - deltaT);
%       phi_oo = phi(t - deltaT - deltaT_o);
%       deltaT_o is the delta time at old time step.
%
%   To find the formula, multiply (I) by (deltaT+deltaT_o), (II) by
%   deltaT and subtracting the first result from the second result.
%
%   Parameters:
% - density = the scalar density.
%
%   Return:
% - LHS = The Left-hand-side d2dt2 contribution to the final linear
% system. 
% - RHS = The Right-hand-side d2dt2 contribution to the final linear
% system.
% 
%   Obs.: 
% - The phi_o, phi_oo stands for phi at old, and old old time respectively.

%% BEG - Initialization
mesh = phi_o.mesh;
LHS = zeros(mesh.numberOfElements,mesh.numberOfElements);
RHS = zeros(mesh.numberOfElements,3);

deltaT = problem.runTime.deltaT();
deltaT_o = problem.runTime.deltaT_o();
a = 2*density / ( (deltaT + deltaT_o) * (deltaT * deltaT_o) );
%% END

%% BEG - Add contribution from time-dependent term.
iElements = 1:mesh.numberOfElements;
idx = sub2ind(size(LHS), iElements, iElements);

LHS(idx) = a * deltaT_o * mesh.volume().data;
RHS(iElements,:) = ...
    bsxfun(@times, a*( ...
          (deltaT + deltaT_o)*phi_o.internalField(iElements).' ...
          - deltaT*phi_oo.internalField(iElements).' ...
      ), mesh.volume().data.');
%% END

end
