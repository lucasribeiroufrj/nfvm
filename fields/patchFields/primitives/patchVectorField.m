classdef patchVectorField
    %patchVectorField An abstraction of a boundary condition.
    %   Detailed explanation goes here
    %
    % TODO:
    % - We should use only a sub-mesh, but should also have a ref. to the
    % main mesh.
    
    properties
        
        vectorField_;
        userName_;
        index_;
        numberOfBFaces_;
        startFace_;
        endFace_;
        mesh_;
        
        conjunctionalCorrection_ = true;
        
        % Stores the discretization coefficient of div(K*gradU). See 
        % TotalLagragian.tractionBoundaryLaplacianOfimpKGradU() for more
        % information.
        divImpKgradU_RHS_;
        divImpKgradU_LHS_;
    end
    
    methods
        
        function obj = patchVectorField(mesh, index, inputVectorField)
            %patchVectorField Construct an instance of 
            %
            %   patchVectorField(mesh, index) is a patch with geometric
            %   boundary information given by
            %   mesh.bondaries(index).
            %
            %   patchVectorField(mesh, index,inputVectorField) is the same 
            %   as above with addition that the field is initialized with
            %   inputVectorField.
            %
            
            boundary = mesh.boundaries(index);
            nFaces = boundary.numberOfBFaces;
            
            % Copy from a boundary.
            if nargin == 2
               
                args{1} = nFaces;
                
            % Initialize the boundary with values from the vectorField.
            elseif nargin == 3 && isa(inputVectorField, 'vectorField')
                
                assert(inputVectorField.numberOfElements == nFaces);
                args{1} = inputVectorField;
            else
                error('No possible construction.')
            end
            
            obj.vectorField_ = vectorField(args{:});
            
            obj.userName_ = boundary.userName;
            obj.index_ = boundary.index;
            obj.numberOfBFaces_ = nFaces;
            obj.startFace_ = boundary.startFace;
            obj.endFace_ = obj.startFace_ + obj.numberOfBFaces_ - 1;
            obj.mesh_ = mesh;
        end
        
        %% BEG - Data member access
        function field = vectorField(obj)
            
            field = obj.vectorField_;
        end
        
        function obj = setVectorField(obj, newValue)
            
            obj.vectorField_ = newValue;
        end
        
        function userName = userName(obj)
            
            userName = obj.userName_;
        end
        
        function index = index(obj)
           
            index = obj.index_;
        end
        
        function nFaces = numberOfFaces(obj)
            
            nFaces = obj.numberOfBFaces_;
        end
        
        function startFace = startFace(obj)
            
            startFace = obj.startFace_;
        end
        
        function endFace = endFace(obj)
           
            endFace = obj.endFace_;
        end
        
        function mesh = mesh(obj)
            
            mesh = obj.mesh_;
        end
        %% END - Data member access
        
        function output = field(obj,i,j)
            
            if nargin == 1
            
                output = obj.vectorField_;
                
            elseif nargin == 2

                output = obj.vectorField_.field(i);

            elseif nargin == 3

                output = obj.vectorField_.field(i,j);

            else
                error('Use one or two parameters');
            end
        end
        
        function newPatch = plus(lhs, rhs)
            
            if isa(lhs, 'patchVectorField') ...
                    && isa(rhs, 'patchVectorField')
                
                newVectorField = ...
                    lhs.vectorField_.plus(rhs.vectorField_);
                
                newPatch = ...
                    patchVectorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newVectorField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = minus(lhs, rhs)
           
            if isa(lhs, 'patchVectorField') ...
                    && isa(rhs, 'patchVectorField')
                
                newVectorField = ...
                    lhs.vectorField_.minus(rhs.vectorField_);
                
                newPatch = ...
                    patchVectorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newVectorField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = mtimes(lhs, rhs)
           
            if length(rhs) == 1 && isfloat(rhs)
                
                newVectorField = lhs.vectorField_.mtimes(rhs);
                
                newPatch = ...
                    patchVectorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newVectorField...
                    );
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newVectorField = rhs.vectorField_.mtimes(lhs);
                
                newPatch = ...
                    patchVectorField...
                    (...
                        rhs.mesh_,...
                        rhs.index_,...
                        newVectorField...
                    );
                
            elseif isa(rhs, 'patchScalarField')
                
                newVectorField = ...
                    lhs.vectorField_.mtimes(rhs.field());
                
                newPatch = ...
                    patchVectorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newVectorField...
                    );
                
            else
                error('Not implemented yet');
            end
        end
        
        function newPatchVectorField = mrdivide(lhs, rhs)
           
            if length(rhs) == 1 && isfloat(rhs)
                
                newVectorField = lhs.vectorField_.mrdivide(rhs);
                
                newPatchVectorField = ...
                    patchVectorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newVectorField...
                    );
                
            elseif isa(rhs, 'patchScalarField')
                
                newVectorField = ...
                    lhs.vectorField_.mrdivide(rhs.field());
                
                newPatchVectorField = ...
                    patchVectorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newVectorField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = mpower(lhs, rhs)
            %mpower Calculate the inner product between two vector fields.
           
            if isa(lhs, 'patchVectorField') && isa(rhs, 'patchVectorField') 
                
                newScalarField = ...
                    lhs.vectorField_.mpower(rhs.vectorField_);
                
                newPatch = ...
                    patchScalarField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newScalarField...
                    );
            else
                ME = MException('patchVectorField:mpower',...
                        'Cannot calculate the inner product.');
                throw(ME);
            end
        end
        
        function newPatch = and(lhs, rhs)
            %and Calculate the dyadic product between two vector fields.
           
            if isa(lhs, 'patchVectorField') && isa(rhs, 'patchVectorField') 
                
                newTensorField = lhs.vectorField_.and(rhs.field());
                
                newPatch = ...
                    patchTensorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
            else
                ME = MException('patchVectorField:and',...
                        'Cannot calculate the dyadic product.');
                throw(ME);
            end
        end
        
        function newPatch = uminus(obj)
           
            field = obj.vectorField_.uminus();
            newPatch = ...
                patchVectorField(obj.mesh_, obj.index_, field);
        end
        
        function newPatch = norm(obj)
           
            normField = obj.vectorField_.norm();
            newPatch = ...
                patchScalarField(obj.mesh_, obj.index_, normField);
        end
    
        function obj = setDivImpKgradUfluxes(obj, LHS, RHS)
            
            if ~isempty(LHS)
                
                obj.divImpKgradU_LHS_ = LHS;
            end
            
            if ~isempty(RHS)
                
                obj.divImpKgradU_RHS_ = RHS;
            end
        end
        
        function [LHS, RHS] = divImpKgradUfluxes(obj)
            
            LHS = obj.divImpKgradU_LHS_;
            RHS = obj.divImpKgradU_RHS_;
        end
        
        function nTotalElements = nTotalElements(obj)
           
            nTotalElements = obj.numberOfBFaces_;
        end
        
        function obj = updateDivImpKgradUfluxes(obj, solidModel)
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            mesh = obj.mesh();
            index = obj.index();
            impK = solidModel.surfaceImpK().boundaryField(index);
            gradU = solidModel.surfaceGradU().boundaryField(index);
            Sf = mesh.Sf().boundaryField(index);
            nf = mesh.nf().boundaryField(index);
            eCN = mesh.eCN().boundaryField(index);
            CNnorm = mesh.CNnorm().boundaryField(index);
            
            [Ef, Tf] = fvmSfDecomposition(Sf, nf, impK, eCN);
            
            gDiff = Ef.norm() / CNnorm;
            FluxVf = gradU * Tf;
            phi_f = obj;
            
            % FluxFf = -FluxCf
            %    RHS = FluxVf + FluxFf*phi_f
            %fluxRHS = FluxVf - FluxCf*phi_f;
            fluxRHS = FluxVf + gDiff*phi_f;
            
            obj.divImpKgradU_LHS_ = -gDiff.field().data;
            obj.divImpKgradU_RHS_ = -(fluxRHS.field().data).';
        end
        
        function obj = correct(obj, volU, volGradU, surfaceGradU) %#ok<INUSD>
            % correct Recompute the field using Tayler series.
            %
            
            owners = obj.mesh().boundaries(obj.index_).owners;
            gradU = tensorField(volGradU.internalField(owners));
            U = vectorField(volU.internalField(owners));
            
            if obj.conjunctionalCorrection_
                
                CN = obj.mesh().CN().boundaryField(obj.index_).field();
                obj.vectorField_ = U + gradU * CN;
            else
                obj.vectorField_ = U;
            end
        end
        
        function obj = setFromRawData(obj, data, fstIdx, sndIdx)
            % setFromRawData(data) return the same field with the new data.
            %   data must accept the indices (1:3, fstIdx:sndIdx) and
            %   fstIdx will be mapped onto the boundary field face of index
            %   1.
            
            obj.vectorField_ = ...
                obj.vectorField_.setFromRawData(data, fstIdx, sndIdx);
        end
    end
end

