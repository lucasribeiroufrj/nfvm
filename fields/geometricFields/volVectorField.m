classdef volVectorField
    %volVectorField Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Interior field
        internalField_;
        
        % Boundary field
        boundaryField_
        
        % FVM Space-time mesh handle 
        mesh_
    end
    
    methods
        
        function obj = volVectorField(mesh, input1, input2)
            %volVectorField Create a geometric field.
            %
            %   volVectorField(mesh) is a geomtric field
            %   with default boundary (calculated/patchVectorField).
            %
            %   volVectorField(mesh, a_vectorField) is the same as
            %   above, except that the field is inialized with values from
            %   a_vectorField.
            %
            %   volVectorField(mesh, boundaryConditions) is a geometric 
            %   field with boundary condition given by the 
            %   boundaryConditions input.
            %
            %   volVectorField(mesh, a_vectorField, a_bVectorField) used to
            %   "mount" another field using the "pieces" a_vectorField and
            %   a_bVectorField. This is the constructor used when defining
            %   overload operator.
            %
            
            fromMeshConstructor = false;
            copyFromRawConstructor = false;
            fromMeshWithBoundaryInfoConstructor = false;
            fromPiecesConstructor = false;
            
            % Contruct from mesh.
            if nargin == 1
                
                args{1} = mesh.numberOfElements;
                fromMeshConstructor = true;
            
            elseif nargin == 2
                
                % Copy constructor from raw data.
                if isa(input1, 'vectorField')
                    
                    nTotalElements = mesh.numberOfTotalElements();
                    assert(input1.numberOfElements == nTotalElements);
                    args{1} = input1.range(1:mesh.numberOfElements());
                    copyFromRawConstructor = true;
                    
                % The second argument must be the boundary conditions set.
                elseif length(input1) == mesh.numberOfBoundaries
                    
                    boundaries = input1;
                    args{1} = mesh.numberOfElements;
                    fromMeshWithBoundaryInfoConstructor = true;
                end
                
            elseif nargin == 3 && isa(input1, 'vectorField') ...
                    && isa(input2, 'bVectorField')
                
                nElements = mesh.numberOfElements();
                assert(input1.numberOfElements() == nElements);
                args{1} = input1;
                fromPiecesConstructor = true;
            else
                error('No possible construction.')
            end
            
            obj.internalField_ = vectorField(args{:});
            obj.mesh_ = mesh;
            
            if fromMeshConstructor == true
                
                obj.boundaryField_ = ...
                    volVectorField.createDefaultBoundary(obj.mesh_);
                
            elseif copyFromRawConstructor
                
                obj.boundaryField_ = ...
                    volVectorField.extractFromVectorField...
                    (...
                        obj.mesh_,...
                        input1...
                    );
                
            elseif fromMeshWithBoundaryInfoConstructor
                
                obj.boundaryField_ = ...
                    volVectorField.extractBoundaryFields...
                    (...
                        obj.mesh_,...
                        boundaries...
                    );
                
            elseif fromPiecesConstructor
                
                obj.boundaryField_ = input2;
            else
                error('No possible construction.')
            end
        end
        
        %% BEG - Data member access
        function mesh = mesh(obj)
            
            mesh = obj.mesh_;
        end
        
        function bVectorField = bVectorField(obj)
            
            bVectorField = obj.boundaryField_;
        end
        
        function output = internalField(obj,i,j)
            
            if nargin == 1

                output = obj.internalField_;
                
            elseif nargin == 2

                output = obj.internalField_.internalField(i);

            elseif nargin == 3

                output = obj.internalField_.internalField(i,j);

            else
                error('Use one or two parameters');
            end
        end
        
        function boundary = boundaryField(obj,b,i,j)
            
            if nargin == 1
                % Hack to support a.boundaryField(1:end)
                boundary = obj.boundaryField_;
            
            elseif nargin == 2
            
                boundary = obj.boundaryField_.boundary(b);
                
            elseif nargin == 3

                boundary = obj.boundaryField_.boundary(b).field(i);

            elseif nargin == 4

                boundary = obj.boundaryField_.boundary(b).field(i,j);
            else
                error('Use one or two parameters');
            end
        end
        
        function nInteriorElements = nInteriorElements(obj)
            % nInteriorElements Returns the # of interior elements
           
            nInteriorElements = obj.internalField_.numberOfElements();
        end
        
        function nTotalElements = nTotalElements(obj)
            % nTotalElements Returns the # of interior + boundary elements
           
            nTotalElements = obj.nInteriorElements() ...
                + obj.boundaryField_.nTotalElements();
        end
        %% END
        
        function value = field(obj)
            % field internal usage function
            %   It is to support, e.g., a.boundaryField(2).field(1:2)
            %
            
            value = obj;
        end
        
        function newVolVectorField = plus(lhs, rhs)
            
            if isa(lhs, 'volVectorField') && isa(rhs, 'volVectorField')
                
                newInternalField = lhs.internalField_ + rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ + rhs.boundaryField_;
                
                newVolVectorField = ...
                    volVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volVectorField:plus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newVolVectorField = minus(lhs, rhs)
            
            if isa(lhs, 'volVectorField') && isa(rhs, 'volVectorField')
                
                newInternalField = lhs.internalField_ - rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ - rhs.boundaryField_;
                
                newVolVectorField = ...
                    volVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volVectorField:minus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newVolVectorField = mtimes(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.internalField_.mtimes(rhs);
                
                newBoundaryField = lhs.boundaryField_ * rhs;
                
                newVolVectorField = ...
                    volVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif length(lhs) == 1 && isfloat(lhs)

                newInternalField = rhs.internalField_.mtimes(lhs);
                
                newBoundaryField = rhs.boundaryField_ * (lhs);
                
                newVolVectorField = ...
                    volVectorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(rhs, 'volScalarField')
                
                newInternalField = ...
                    lhs.internalField_.mtimes(rhs.internalField());
                
                newBoundaryField = lhs.boundaryField_ * rhs.boundaryField();
                
                newVolVectorField = ...
                    volVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volVectorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newVolVectorField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 ...
                    && (isfloat(rhs) || isa(rhs, 'volScalarField') )
                
                newInternalField = ...
                    lhs.internalField_.mrdivide(rhs.internalField());
                
                newBoundaryField = lhs.boundaryField_ / rhs.boundaryField();
                
                newVolVectorField = ...
                    volVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volVectorField:mrdivide',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newVolScalarField = mpower(lhs, rhs)
            
            if isa(lhs, 'volVectorField') && isa(rhs, 'volVectorField')
                
                newInternalField = ...
                    lhs.internalField_.mpower(rhs.internalField_);
                
                newBoundaryField = lhs.boundaryField_ ^ rhs.boundaryField_;
                
                newVolScalarField = ...
                    volScalarField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volVectorField:mpower',...
                        'Cannot calculate the inner product.');
                throw(ME);
            end
        end
        
        function newVolTensorField = and(lhs, rhs)
            %and Tensorial product between two vectors.
            
            if isa(lhs, 'volVectorField') ...
                    && isa(rhs, 'volVectorField')
                
                newInternalField = lhs.internalField_ & rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ & rhs.boundaryField_;
                
                newVolTensorField = ...
                    volTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volVectorField:and',...
                        'Cannot calculate the dyadic product.');
                throw(ME);
            end
        end
        
        function newVectorField = uminus(obj)
           
            newInternalField = obj.internalField_.uminus();
                
            newBoundaryField = obj.boundaryField_.uminus();

            newVectorField = ...
                volVectorField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newVolScalarField = norm(obj)
           
            newInternalField = obj.internalField_.norm();
                
            newBoundaryField = obj.boundaryField_.norm();

            newVolScalarField = ...
                volScalarField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function obj = updateDivImpKgradUfluxes(obj, solidModel)
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            obj.boundaryField_ = ...
                obj.boundaryField_.updateDivImpKgradUfluxes(solidModel);
        end
        
        function obj = ...
                correctBoundaryConditions...
                (...
                    obj,...
                    volU,...
                    volGradU,...
                    surfaceGradU...
                )
            
            obj.boundaryField_ = ...
                obj.boundaryField_.correctBoundaryConditions...
                (...
                    volU,...
                    volGradU,...
                    surfaceGradU...
                );
        end
        
        function obj = setFromRawData(obj, data, fstIdx, sndIdx)
            % setFromRawData(data) return the same field with the new data.
            %   data must accept the indices (1:3, startIdx:endIdx) and
            %   fstIdx will be mapped onto the cell of index 1.
            
            nElements = obj.nTotalElements();
            
            if fstIdx > sndIdx || fstIdx < 1 || ...
                    (sndIdx - fstIdx + 1) > nElements
                error('invalid indexing');
            end
            
            rgtIdx = obj.nInteriorElements();
            
            if fstIdx <= rgtIdx
                
                if sndIdx <= rgtIdx
                    
                    k = sndIdx;
                    left = 0;
                else
                    left = sndIdx - rgtIdx;
                    k = rgtIdx;
                end
                
                obj.internalField_ = ...
                    obj.internalField_.setFromRawData(data, fstIdx, k);
                fstIdx = k + 1;
            else
                left = sndIdx - rgtIdx;
            end
            
            if left > 0
                
                obj.boundaryField_ = ...
                    obj.boundaryField_.setFromRawData(data,fstIdx,sndIdx); 
            end
        end
        
        function obj = setInternalData(obj, newValue)
            
            obj.internalField_ = obj.internalField_.setData(newValue);
        end
    end
    
    methods (Static)
        
        function bfield = extractFromVectorField(mesh, inputVecField)
            
            nBoundaries = mesh.numberOfBoundaries;
            bfield = bVectorField(nBoundaries);
            
            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                startIndex = ...
                    boundary.startFace - mesh.numberOfInteriorFaces...
                    + mesh.numberOfElements;
                nFaces = boundary.numberOfBFaces;
                endIndex = startIndex + nFaces - 1;
                segment = inputVecField.subField(startIndex:endIndex);
                
                pField = patchVectorField(mesh, boundary.index, segment);
                bfield = bfield.setPatch(b,pField);
            end
        end
        
        function bField = createDefaultBoundary(mesh)
            
            nBoundaries = mesh.numberOfBoundaries;
            bField = bVectorField(nBoundaries);

            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                pField = patchVectorField(mesh, boundary.index);
                bField = bField.setPatch(b,pField);
            end
        end
        
        function bField = extractBoundaryFields(mesh, boundaries) %#ok<INUSL>

            nBoundaries = length(boundaries);
            bField = bVectorField(nBoundaries);
            
            for b = 1:nBoundaries
                
                boundary = boundaries{b};
                kind = boundary.kind;
                radical = 'PatchVectorField';
                pField = eval([kind, radical, '(mesh, boundary)']);
                bField = bField.setPatch(b,pField);
            end
        end
    end
end

