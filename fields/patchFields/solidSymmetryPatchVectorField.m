classdef solidSymmetryPatchVectorField < patchVectorField
    %solidSymmetryPatchVectorField An abstraction of a boundary condition.
    %   % The symmetry boundary condition [eq.3.46, Maneeratana - 2000;
    %   eq.3.45, Muzaferija - 1994]. This scheme implementats a new 
    %   correction for non-conjuntional meshes.
    %
    % TODO: 
    %  - To be compatible with solids4Foam:
    % NAME
    %{
    %    type            solidSymmetry;
    %    patchType       symmetryPlane;
    %    value           uniform (0 0 0);
    %}
    
    properties
        
        Ef_;
        Tf_;
        
        d_fprime_f_;
        d_C_Cprime_;
        d_Cp_C2pNorm_;
    end
    
    methods
        
        function obj = ...
                solidSymmetryPatchVectorField...
                (...
                    mesh,...
                    boundary,...
                    inputVectorField...
                )
            %solidSymmetryPatchVectorField Construct an instance of this
            % class Detailed explanation goes here
            
            args{1} = mesh;
            args{2} = boundary.index;
            
            if nargin == 3
                
                args{3} = inputVectorField;
            end
            
            obj@patchVectorField(args{:});
            
            obj = obj.initialized();
        end
        
        function obj = initialized(obj)
            %initialized Creates some cache.
            
            %% Needed for correct() and updateDivImpKgradUfluxes() methods.
            mesh = obj.mesh();
            index = obj.index();
            owners = mesh.boundaries(index).owners;
            r_C = vectorField(mesh.r_C().internalField(owners));
            r_Cf = mesh.r_C().boundaryField(index).field();
            d_CN = mesh.CN().boundaryField(index).field();
            nf = mesh.nf().boundaryField(index).field();
            
            r_Cfprime = (d_CN ^ nf)*nf;
            obj.d_fprime_f_ = (r_C + r_Cfprime) - r_Cf;
            %%
            
            %% BEG - Needed for updateDivImpKgradUfluxes.
            d_CN_perp = r_Cfprime; % To correspond to the manual.
            obj.d_C_Cprime_ = d_CN - d_CN_perp;             % (1)
            %phi_Cprime = phi_C + gradPhi_C*d_C_Cprime_;    % (2)
            r_Cprime = r_C + obj.d_C_Cprime_;               % (3)
            r_C2prime = r_Cprime + 2*d_CN_perp;             % (4)
            %phi_C2p = phi_Cprime - 2*(phi_Cprime.'*nf)*nf  % (5)
            d_Cp_C2p = r_C2prime - r_Cprime;                % (6)
            
            Sf = mesh.Sf().boundaryField(index).field();
            Sfnorm = mesh.SfNorm().boundaryField(index).field();
            
            d_Cp_C2pNorm = d_Cp_C2p.norm();
            e_Cp_C2p = d_Cp_C2p / d_Cp_C2pNorm;
            [Ef, Tf] = fvmSfDecompositionNoImpk(Sf,Sfnorm,nf,e_Cp_C2p);
            obj.Ef_ = Ef;
            obj.Tf_ = Tf;
            obj.d_Cp_C2pNorm_ = d_Cp_C2pNorm;
            %% END
        end
        
        function obj = updateDivImpKgradUfluxes(obj, solidModel)
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            mesh = obj.mesh();
            index = obj.index();
            nf = mesh.nf().boundaryField(index).field();
            owners = mesh.boundaries(obj.index_).owners;
            U_C = vectorField(solidModel.volU().internalField(owners));
            
            gradU_C = ...
                tensorField(solidModel.volGradU().internalField(owners));
            
            d_C_Cprime = obj.d_C_Cprime_;
            
            U_Cprime = U_C + gradU_C*d_C_Cprime;    % (2)
            proj = 2*(U_Cprime ^ nf);
            U_C2prime = U_Cprime - proj*nf;         % (5)
            
            gradU_f = ...
                solidModel.surfaceGradU().boundaryField(index).field();
            
            impK = solidModel.surfaceImpK().boundaryField(index).field();
            Tfprime = impK.' * obj.Tf_;
            Efprime = impK.' * obj.Ef_;
            FluxVf = gradU_f * Tfprime;
            
            gDiff = norm(Efprime) / obj.d_Cp_C2pNorm_;
            
            FluxFf = gDiff;
            fluxRHS = FluxFf*U_C2prime + FluxVf;

            divImpKgradU_LHS_ = -FluxFf.data;
            divImpKgradU_RHS_ = -(fluxRHS.data).';
            
            obj = ...
                obj.setDivImpKgradUfluxes...
                ( ...
                    divImpKgradU_LHS_,...
                    divImpKgradU_RHS_...
                );
        end
        
        % ATTENTION to the inclusion of the surfaceGradU
        function obj = ...
                correct(obj, volU, volGradU, surfaceGradU) %#ok<INUSL>
            %correct Re-calculates this boundary field.
            
            mesh = obj.mesh();
            index = obj.index();
            nf = mesh.nf().boundaryField(index).field();
            owners = mesh.boundaries(obj.index_).owners;
            
            U_C = vectorField(volU.internalField(owners));
            U_f = U_C - (U_C ^ nf)*nf;
            
            if obj.conjunctionalCorrection_
                
                gradU_f = surfaceGradU().boundaryField(index).field();
                U_f = U_f - gradU_f*obj.d_fprime_f_;
            end
            
            obj = obj.setVectorField(U_f);
        end
    end
end

