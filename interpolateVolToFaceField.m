function facePhi = ...
    interpolateVolToFaceField(phi, method, input1, input2, input3)
%interpolateVolToFaceField Given a element field it returns the
% interpolated face field.
%   

% FEATURES:
% - The solidSymmetric boundary conditions scheme adopted implements
% a correction for non-conjunctional meshes.
% - The solidTraction boundary condition is implemented using 
% options 2: 1) and 3) in [Pag. 278 - Moukalled] for internal faces and 
% using taylor to calculate the value of phi at the face centroid.

% ATTENTION:
% - Only two boundary conditions are implemented, see above.

conjunctional_correction = true;
isOptimized = true;

mesh = phi.mesh();

if strcmp(method, 'Compact-Stencil')
    
    if nargin == 4
        
        volGradPhi = input1;
        surfaceGradPhi = input2;
        
        %% Calculate phi at internal face centroids using compact
        % stencil method (option 2: 1.) [Pag. 278 - Moukalled].
        
        %%
        facePhi = surfaceVectorField(mesh);
        indices = 1:mesh.numberOfInteriorFaces;
        owners = [mesh.faces(indices).iOwner];
        neighs = [mesh.faces(indices).iNeighbour];
        facePhi.internalField(1:end) = ...
            (phi.internalField(owners) + phi.internalField(neighs))/2;
        facePhi = ...
            facePhi.correctBoundaryConditions...
            (...
                phi,...
                volGradPhi,...
                surfaceGradPhi...
            );
        return;
        %%

        %% Calculate phi at the boundary face centroids.
        for iBoundary=1:mesh.numberOfBoundaries

            boundary = mesh.boundaries(iBoundary);
            
            boundary.kind = 'solidTraction';

            if strcmp(boundary.type,'empty') == true
                continue;
            end
            
            bName = boundary.userName;

            startFace = boundary.startFace;
            endFace = startFace + boundary.numberOfBFaces - 1;
            
            if strcmp(boundary.kind, 'Dirichlet')
                
                for iFace=startFace:endFace

                    phi_f = BOUNDARY.(bName)(iFace - startFace + 1);

                    facePhi(iFace) = phi_f;
                end

            elseif strcmp(boundary.kind, 'Neumann') ...
                    || strcmp(boundary.kind, 'ZeroFluxNeumann')
                
                for iFace=startFace:endFace
                    
                    iOwner = mesh.faces(iFace).iOwner;
                    facePhi(iFace) = phi(iOwner);
                end
                
            elseif strcmp(boundary.kind, 'solidTraction')
                
                if ~isOptimized
                    for iFace=1:boundary.numberOfBFaces

                        face = mesh.faces(startFace + iFace - 1);
                        iOwner = face.iOwner;

                        if conjunctional_correction

                            gradU = volGradPhi.internalField(iOwner);
                            facePhi.boundaryField(iBoundary).field(iFace) = ...
                                phi.internalField(iOwner) + gradU*face.CN;
                        else
                            facePhi.boundaryField(iBoundary).field(iFace) = ...
                                phi.internalField(iOwner);
                        end
                    end
                else
                    endFace = startFace + boundary.numberOfBFaces - 1;
                    owners = [mesh.faces(startFace:endFace).iOwner];
                    
                    gradU = tensorField(volGradPhi.internalField(owners));
                    Phi = vectorField(phi.internalField(owners));
                    
                    if conjunctional_correction
                        
                        CN = vectorField([mesh.faces(startFace:endFace).CN]);
                        newFacePhi = Phi + gradU * CN;
                    else
                        newFacePhi = Phi;
                    end
                    
                    facePhi.boundaryField(iBoundary).field(1:end) = newFacePhi(1:end);
                end
                
            elseif strcmp(boundary.kind, 'solidSymmetry')
                
                % ATTENTION: no non-conformal correction for solidSymmetric boundary
                for iFace=startFace:endFace
                    
                    face = mesh.faces(iFace);
                    iOwner = face.iOwner;
                    n = face.Sf / norm(face.Sf);
                    phi_C = phi(iOwner);
                    facePhi(iFace) = phi_C ...
                        - (phi_C.'*n)*n;
                    
                    %% Non-conformal correction
                    if conjunctional_correction
                        r_C = mesh.elements(iOwner).centroid;
                        r_Cfprime = (face.CN.'*n)*n;
                        d_fprime_f = (r_C + r_Cfprime) - face.centroid;
                        facePhi(iFace) = facePhi(iFace) ...
                            - FIELDS.surfacePhiGrad(iFace)*d_fprime_f;
                    end
                end
            else
                ME = MException(...
                    'Error: Unknown boundary condition type: %s', ...
                    boundary.kind);
                throw(ME)
            end
        end
    else
        
        facePhi = input1;
        volGradPhi = input2;
        
        %% BEG - Calculate phi at internal face centroids using compact
        % stencil method (option 2: 3) [Pag. 278 - Moukalled].
        for iFace=1:mesh.numberOfInteriorFaces
            
            face = mesh.faces(iFace);
            iOwner = face.iOwner;
            iNeigh = face.iNeighbour;
            phi_f_prime = facePhi.internalField(iFace);
            grad_C = volGradPhi.internalField(iOwner);
            grad_F = volGradPhi.internalField(iNeigh);
            rf = face.centroid;
            rC = mesh.elements(iOwner).centroid;
            rF = mesh.elements(iNeigh).centroid;

            facePhi.internalField(iFace) = phi_f_prime ...
                + 0.5*(grad_C + grad_F)*(rf - 0.5*(rC + rF));
        end
        %%

        %% Calculate phi at the boundary face centroids.
        for iBoundary=1:mesh.numberOfBoundaries
            
            boundary = mesh.boundaries(iBoundary);

            boundary.kind = 'solidTraction';
            
            if strcmp(boundary.type,'empty') == true
                continue;
            end

            startFace = boundary.startFace;
            endFace = startFace + boundary.numberOfBFaces - 1;
                
            if strcmp(boundary.kind, 'Dirichlet')
                
                % No need, since it was already computed during the first
                % iteration
            
            elseif strcmp(boundary.kind, 'solidTraction')
                
                if ~isOptimized
                    for iFace=1:boundary.numberOfBFaces

                        face = mesh.faces(startFace + iFace - 1);
                        iOwner = face.iOwner;

                        if conjunctional_correction

                            gradU = volGradPhi.internalField(iOwner);
                            facePhi.boundaryField(iBoundary).field(iFace) = ...
                                phi.internalField(iOwner) + gradU*face.CN;
                        else
                            facePhi.boundaryField(iBoundary).field(iFace) = ...
                                phi.internalField(iOwner);
                        end
                    end
                else
                    endFace = startFace + boundary.numberOfBFaces - 1;
                    owners = [mesh.faces(startFace:endFace).iOwner];
                    
                    gradU = tensorField(volGradPhi.internalField(owners));
                    Phi = vectorField(phi.internalField(owners));
                    
                    if conjunctional_correction
                        
                        CN = vectorField([mesh.faces(startFace:endFace).CN]);
                        newFacePhi = Phi + gradU * CN;
                    else
                        newFacePhi = Phi;
                    end
                    
                    facePhi.boundaryField(iBoundary).field(1:end) = newFacePhi(1:end);
                end
                
            elseif strcmp(boundary.kind, 'solidSymmetry')
                
                % No need, since it was already computed during the first
                % iteration
            else
                ME = MException('Error: Unknown boundary type %s', ...
                    boundary.kind);
                throw(ME)
            end
        end
    end
else
    ME = MException('Error: There is no interpolation called %s',method);
    throw(ME)
end

end

