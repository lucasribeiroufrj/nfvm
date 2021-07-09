function volGradPhi = ...
    fvc_grad(phi, Sf, config, volGradPhi_o, surfaceGradPhi_o)
%fvc_grad Summary of this function goes here
%   Detailed explanation goes here

isToCorrectBoundary = false;

if nargin == 5
    isToCorrectBoundary = true;
end

if strcmp(config.gradInterpolationScheme, 'Gauss-Green')
% Calculate gradient at cell centroids using 
% compact stencil method (option 2) [Pag. 278 - Moukalled]
    
    for gradIterComputation=1:config.maxGradIterComputation

        if gradIterComputation == 1

            mesh = phi.mesh();
            facePhi = surfaceVectorField...
                (...
                    mesh,...
                    vectorField(mesh.numberOfInteriorFaces),...
                    phi.boundaryField()...
                );
            indices = 1:mesh.numberOfInteriorFaces();
            owners = [mesh.faces(indices).iOwner];
            neighs = [mesh.faces(indices).iNeighbour];
            
            facePhi = ...
                facePhi.setInternalData...
                (...
                    (phi.internalField(owners) ...
                    + phi.internalField(neighs))/2 ...
                );
            
            if isToCorrectBoundary
                
                facePhi = ...
                    facePhi.correctBoundaryConditions...
                    (...
                        phi,...
                        volGradPhi_o,...
                        surfaceGradPhi_o...
                    );
            end
        else
            
            mesh = phi.mesh();
            
            phi_f_prime = facePhi.internalField();
            
            iOwners = mesh.iOwners();
            iNeighs = mesh.iNeighs();
            
            grad_C = tensorField(volGradPhi.internalField(iOwners));
            grad_F = tensorField(volGradPhi.internalField(iNeighs));
            
            rf = mesh.r_Cf.internalField();
            rC = vectorField(mesh.r_C_.internalField(iOwners));
            rF = vectorField(mesh.r_C_.internalField(iNeighs));
            
            facePhi = ...
                facePhi.setInternalField...
                (...
                    phi_f_prime + ...
                    0.5*(grad_C + grad_F)*(rf - 0.5*(rC + rF))...
                );
            
            if isToCorrectBoundary
                
                facePhi = ...
                    facePhi.correctBoundaryConditions...
                    (...
                        phi,...
                        volGradPhi,...
                        surfaceGradPhi_o...
                    );
            end
        end

        volGradPhi = fvc_GaussGreenGrad(facePhi, Sf);
    end
    
elseif strcmp(config.gradInterpolationScheme, 'Least-Square')
    
    volGradPhi = fvc_LeastSquareGrad(phi);
    
else
    fprintf('There no such gradInterpolationScheme: %s\n', ...
        config.gradInterpolationScheme);
    exit
end

end

