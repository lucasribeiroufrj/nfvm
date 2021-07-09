function [LHS, RHS] = laplacianAtBoundary(...
    impK, surfaceGradPhi, boundary, phi)
%laplacianAtBoundary Summary of this function goes here
%   Detailed explanation goes here

global FIELDS
global BOUNDARY

NEW_SYMMETRY_SCHEME = true;
USE_MY_SYMMETRIC_SOLUTION = false;

if USE_MY_SYMMETRIC_SOLUTION
    conjunctional_correction = true;
end

mesh = surfaceGradPhi.mesh();

LHS = zeros(mesh.numberOfElements,mesh.numberOfElements);
RHS = zeros(mesh.numberOfElements,3);

Sf = [mesh.faces(:).Sf];

% Prescibed displacement
% TODO - rigidSphereCollision should be a Dirichlet class!
if strcmp(boundary.kind, 'Dirichlet') ...
        || strcmp(boundary.kind, 'rigidSphereCollision') ...
        || strcmp(boundary.kind, 'zeroDisplacement')
    
    startFace = boundary.startFace;
    endFace = startFace + boundary.numberOfBFaces - 1;
    bName = boundary.userName;

    for iFace=startFace:endFace
        
        face = mesh.faces(iFace);
        [Ef, Tf] = fvmSfDecomposition(mesh,Sf,iFace,impK);
        gDiff = norm(Ef) / norm(face.CN);
        FluxCf = -gDiff;
        FluxFf = -FluxCf;
        gradPhi = surfaceGradPhi(iFace);
        FluxVf = gradPhi * Tf;                
        iOwner = face.iOwner;

        LHS(iOwner,iOwner) = LHS(iOwner,iOwner) + FluxCf;

        phi_f = BOUNDARY.(bName)(:, iFace - startFace + 1);

        RHS(iOwner,:) = RHS(iOwner,:) - FluxVf.' - FluxFf*phi_f.';
    end
    
elseif strcmp(boundary.kind, 'solidSymmetry')

    % This discretization works with skewed and non-conformal meshes. It 
    % also supports general diffusion term.
    startFace = boundary.startFace;
    endFace = startFace + boundary.numberOfBFaces - 1;

    for iFace=startFace:endFace
        
        % Perp stands for perpendicular
        face = mesh.faces(iFace);
        
        if NEW_SYMMETRY_SCHEME
            
            d_CN = face.CN;
            n = face.Sf/norm(face.Sf);
            d_CN_perp = (d_CN.'*n)*n;
            iOwner = face.iOwner;
            element = mesh.elements(iOwner);
            r_C = element.centroid;
            phi_C = phi(iOwner);
            gradPhi_C = FIELDS.volPhiGrad(iOwner);
            d_C_Cprime = d_CN - d_CN_perp;                  % (1)
            phi_Cprime = phi_C + gradPhi_C*d_C_Cprime;      % (2)
            r_Cprime = r_C + d_C_Cprime;                    % (3)
            r_C2prime = r_Cprime + 2*d_CN_perp;             % (4)
            phi_C2p = phi_Cprime - 2*(phi_Cprime.'*n)*n;    % (5)
            d_Cp_C2p = r_C2prime - r_Cprime;                % (6)
            
            [Ef, Tf] = fvmSfDecompositionForGhost(Sf, iFace, d_Cp_C2p);
            gDiff = norm(Ef) / norm(d_Cp_C2p);
            FluxCf = -gDiff;
            FluxFf = -FluxCf;
            gradPhi = surfaceGradPhi(iFace);
            FluxVf = gradPhi * Tf;
            phi_b = phi_C2p;
        else
        
            % Perp stands for perpendicular
            [Ef, Tf] = fvmSfDecomposition(mesh,Sf,iFace,impK);
            n = face.Sf/norm(face.Sf);
            gDiff = norm(Ef) / norm(face.CN);
            FluxCf = -gDiff;
            FluxFf = -FluxCf;
            gradPhi = surfaceGradPhi(iFace);
            FluxVf = gradPhi * Tf;

            if USE_MY_SYMMETRIC_SOLUTION
                % Check [eq.3.46, Maneeratana - 2000; Muzaferija - 1994].
                iOwner = face.iOwner;
                phi_C = phi(iOwner);
                phi_b = phi_C - (phi_C.'*n)*n;

                % Non-conformal correction. We just use taylor from f to f' 
                % (to understand this notation f' check [fig. 9.3c,
                % Moukalled - 2015]).
                if conjunctional_correction
                    CN_perp = (face.CN.'*n)*n;
                    d_f_fprime = CN_perp - face.CN;
                    phi_b = phi_b - gradPhi*d_f_fprime;
                end
            else
                % BEG - Non-conformal correction.
                iOwner = face.iOwner;
                phi_C = phi(iOwner);
                gradPhi_C = FIELDS.volPhiGrad(iOwner);
                d_Cprime_C = face.CN - (face.CN.'*n)*n;
                phi_Cprime = phi_C + gradPhi_C*d_Cprime_C;
                % END

                % Check [eq.3.46, Maneeratana - 2000; Muzaferija - 1994].
                phi_b = phi_Cprime - (phi_Cprime.'*n)*n;
            end
            
        end
            
        LHS(iOwner, iOwner) = LHS(iOwner, iOwner) + FluxCf;
        RHS(iOwner,:) = RHS(iOwner,:) - FluxFf*phi_b.' - FluxVf.';
    end
    
elseif strcmp(boundary.kind, 'solidTraction') % Neumann
    
    startFace = boundary.startFace;
    endFace = startFace + boundary.numberOfBFaces - 1;
    bName = boundary.userName;
    
    for iFace=startFace:endFace
        
        face = mesh.faces(iFace);
        S = face.Sf;
        idx = iFace - startFace + 1;
        gradPhi = surfaceGradPhi.boundaryField(idx).field(idx);
        resultantforceAtBoundary = 1;%BOUNDARY.(bName)(iFace - startFace + 1);
        flux = resultantforceAtBoundary + impK.internalField(iFace)*gradPhi*S;
        iOwner = face.iOwner;
        RHS(iOwner,:) = RHS(iOwner,:) - flux';
    end

elseif strcmp(boundary.kind, 'tractionFree') % Neumann
    
    startFace = boundary.startFace;
    endFace = startFace + boundary.numberOfBFaces - 1;
    
    new = true;
    if new
        A = superclasses(class(impK)); 
        className = A{:};
        Sp = vectorField([ mesh.faces(startFace:endFace).Sf ]);
        gradPhi_p = tensorField(surfaceGradPhi(startFace:endFace));
        impKp = eval([className,'(impK(startFace:endFace))']);
        flux = gradPhi_p*(impKp.'*Sp);
        RHS([mesh.faces(startFace:endFace).iOwner],:) = ...
            RHS([mesh.faces(startFace:endFace).iOwner],:) + flux.data';
    else
        for iFace=startFace:endFace

            face = mesh.faces(iFace);
            S = face.Sf;
            gradPhi = surfaceGradPhi(iFace);
            flux = impK(iFace)*gradPhi*S;
            iOwner = face.iOwner;
            RHS(iOwner,:) = RHS(iOwner,:) - flux';
        end
    end
else
    ME = MException('Error: Unknown boundary condition type: %s', ...
        boundary.kind);
    throw(ME)
end

end
