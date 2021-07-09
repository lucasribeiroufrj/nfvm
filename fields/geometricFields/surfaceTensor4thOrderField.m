classdef surfaceTensor4thOrderField < tensor4thOrderField
    %surfaceTensor4thOrderField Summary of this class goes here
    %   Detailed explanation goes here
    %
    % TODO:
    % - We should use only a sub-mesh, but should also have a ref. to the
    % main mesh.
    
    properties
        mesh_
        
        % Boundary field
        bTensor4thOrderField_
        
        % Number of interior elements
        nInteriorElements_
        
        % Number of interior + boundary elements
        nTotalElements_
    end
    
    methods
        function obj = surfaceTensor4thOrderField(mesh, input1, input2)
            %surfaceTensor4thOrderField Create a geometric field.
            %
            %   surfaceTensor4thOrderField(mesh) is a geomtric field
            %   with default boundary (calculated/patchTensor4thOrderField).
            %
            %   surfaceTensor4thOrderField(mesh, a_tensor4thOrderField)
            %   is the same as above, except that the field is inialized
            %   with values from a_tensor4thOrderField.
            %
            %   surfaceTensor4thOrderField(mesh, boundaries) is geometric
            %   field with boundary information given by the boundaries
            %   array.
            %
            %   surfaceTensor4thOrderField
            %   (...
            %       mesh,...
            %       a_tensor4thOrderField,...
            %       a_bTensor4thOrderField...
            %   ) 
            %   used to "mount" another field using the "pieces" 
            %   a_tensor4thOrderField and a_bTensor4thOrderField. This is
            %   the constructor used when defining overload operator.
            %
            
            fromMeshConstructor = false;
            copyFromTensorFieldConstructor = false;
            fromMeshWithBoundaryInfoConstructor = false;
            fromPiecesConstructor = false;
            
            % Contruct from mesh.
            if nargin == 1
                
                args{1} = mesh.numberOfInteriorFaces;
                fromMeshConstructor = true;
            
            % Copy constructor from raw data.
            elseif nargin == 2
                
                if isa(input1, 'tensor4thOrderField')
                
                    nElements = mesh.numberOfInteriorFaces;
                    internalTensorField = ...
                        tensor4thOrderField(input1(1:nElements));
                    args{1} = internalTensorField;
                    copyFromTensorFieldConstructor = true;
                    
                % The second argument must be the boundary conditions set.
                elseif isfield(input1, 'index') ...
                        && isfield(input1, 'kind') ...
                        && length(input1) == mesh.numberOfBoundaries
                    
                    boundaries = input1;
                    args{1} = mesh.numberOfInteriorFaces;
                    fromMeshWithBoundaryInfoConstructor = true;
                else
                    error('No possible construction.')
                end
                
            elseif nargin == 3 && isa(input1, 'tensor4thOrderField') ...
                    && isa(input2, 'bTensor4thOrderField')
                
                args{1} = input1;
                fromPiecesConstructor = true;
            else
                error('No possible construction.')
            end
            
            obj@tensor4thOrderField(args{:});
            obj.mesh_ = mesh;
            
            if fromMeshConstructor == true
                
                obj.bTensor4thOrderField_ = ...
                    surfaceTensor4thOrderField ...
                    .createDefaultBoundary(obj.mesh_);
                
            elseif copyFromTensorFieldConstructor
                
                obj.bTensor4thOrderField_ = ...
                    surfaceTensor4thOrderField.extractFromTensorField...
                    (...
                        obj.mesh_,...
                        input1...
                    );
                
            elseif fromMeshWithBoundaryInfoConstructor
                
                obj.bTensor4thOrderField_ = ...
                    surfaceTensor4thOrderField.extractBoundaryFieds...
                    (...
                        obj.mesh_,...
                        boundaries...
                    );
                
            elseif fromPiecesConstructor
                    
                obj.bTensor4thOrderField_ = input2;
            else
                error('No possible construction.')
            end
            
            obj.nInteriorElements_ = obj.mesh_.numberOfInteriorFaces;
            obj.nTotalElements_ = obj.mesh_.numberOfFaces;
        end
        
        function output = internalField(obj,i,j)
            
            if nargin == 1
                % Hack to support a.internalField(1:end)
                output = tensor4thOrderField(obj);
                
            elseif nargin == 2

                output = internalField@tensor4thOrderField(obj,i);

            elseif nargin == 3

                output = internalField@tensor4thOrderField(obj,i,j);

            else
                error('Use one or two parameters');
            end
        end
        
        function boundary = boundaryField(obj,b,i,j)
            
            if nargin == 1

                boundary = obj.bTensor4thOrderField_;
            
            elseif nargin == 2
            
                boundary = obj.bTensor4thOrderField_.boundary(b);
                
            elseif nargin == 3

                boundary = obj.bTensor4thOrderField_.boundary(b).field(i);

            elseif nargin == 4

                boundary = obj.bTensor4thOrderField_.boundary(b).field(i,j);
            else
                error('Use one or two parameters');
            end
        end
        
        function value = field(obj)
            % field Internal usage function
            %   It is to support, e.g., a.boundaryField(2).field(1:2)
            %
            
            value = obj;
        end
        
        function newSurfaceTensor4thOrderField = plus(lhs, rhs)
            
            if isa(lhs, 'surfaceTensor4thOrderField') ...
                    && isa(rhs, 'surfaceTensor4thOrderField')
                
                newInternalField = lhs.plus@tensor4thOrderField(rhs);
                
                newBoundaryField = lhs.bTensor4thOrderField_...
                   + rhs.bTensor4thOrderField_;
                
                newSurfaceTensor4thOrderField = ...
                    surfaceTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                errorOperation('Cannot sum operands.');
            end
        end
        
        function newSurfaceTensor4thOrderField = minus(lhs, rhs)
            
            if isa(lhs, 'surfaceTensor4thOrderField') ...
                    && isa(rhs, 'surfaceTensor4thOrderField')
                
                newInternalField = lhs.minus@tensor4thOrderField(rhs);
                
                newBoundaryField = lhs.bTensor4thOrderField_ ...
                    - rhs.bTensor4thOrderField_;
                
                newSurfaceTensor4thOrderField = ...
                    surfaceTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else                    
                errorOperation('Cannot subtract operands.');
            end
        end
        
        function newGeometricField = mtimes(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.mtimes@tensor4thOrderField(rhs);
                
                newBoundaryField = lhs.bTensor4thOrderField_ * rhs;
                
                newGeometricField = ...
                    surfaceTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newInternalField = rhs.mtimes@tensor4thOrderField(lhs);
                
                newBoundaryField = rhs.bTensor4thOrderField_ * lhs;
                
                newGeometricField = ...
                    surfaceTensor4thOrderField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                errorOperation('Cannot multiply operands.');
            end
        end
        
        function newSurfaceTensor4thOrderField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.mrdivide@tensor4thOrderField(rhs);
                
                newBoundaryField = lhs.bTensor4thOrderField_ / rhs;
                
                newSurfaceTensor4thOrderField = ...
                    surfaceTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                errorOperation('Cannot divide operands.');
            end
        end
        
        function newSurfaceField = doubleDot(lhs, rhs)
            %doubleDot double constraction like lhs:rhs.
            %   doubleDot(lhs, rhs) returns a surfaceTensorField. This
            %   tensor field being the result of the double contraction:
            %       lhs_{abcd}*rhs_{cd}.
            %   lhs: any surfaceTensor4thOrderField.
            %   rhs: any surfaceTensorField.
            
            if isa(lhs, 'surfaceTensor4thOrderField') ...
                    && isa(rhs, 'surfaceTensorField')
                
                newInternalField = ...
                    lhs.doubleDot@tensor4thOrderField(rhs);
                
                newBoundaryField = lhs.bTensor4thOrderField_...
                    .doubleDot(rhs.boundaryField());
                
                newSurfaceField = ...
                    surfaceTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                errorOperation('Cannot double dot operands.');
            end
        end
        
        function mesh = mesh(obj)
            
            mesh = obj.mesh_;
        end
        
        function bTensor4thOrderField = bTensor4thOrderField(obj)
            
            bTensor4thOrderField = obj.bTensor4thOrderField_;
        end
        
        function nInteriorElements = nInteriorElements(obj)
           
            nInteriorElements = obj.nInteriorElements_;
        end
        
        function nTotalElements = nTotalElements(obj)
           
            nTotalElements = obj.nTotalElements_;
        end
        
        function Identity = identity(obj)
            
            Identity = obj.mesh_.surfaceI();
        end
        
        function obj = copyBoundaryFrom(obj, other)
            %copyBoundaryFrom Assumes the boundaies are compatible.
            
            obj.bTensor4thOrderField_ = other.bTensor4thOrderField_;
        end
    end
    
    methods(Static)
       
        function bfield = extractFromTensorField(mesh, inputTenField)
            
            nBoundaries = mesh.numberOfBoundaries;
            bfield = bTensor4thOrderField(nBoundaries);
            
            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                startIndex = boundary.startFace;
                nFaces = boundary.numberOfBFaces;
                endIndex = startIndex + nFaces - 1;
                segment = ...
                    tensor4thOrderField(inputTenField(startIndex:endIndex));
                
                pField = ...
                    patchTensor4thOrderField(mesh, boundary.index, segment);
                bfield = bfield.setPatch(b,pField);
            end
        end
        
        function bField = createDefaultBoundary(mesh)
            
            nBoundaries = mesh.numberOfBoundaries;
            bField = bTensor4thOrderField(nBoundaries);

            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                pField = patchTensor4thOrderField(mesh, boundary.index);
                bField = bField.setPatch(b,pField);
            end
        end
        
        function bField = extractBoundaryFields(mesh, boundaries) %#ok<INUSL>

            nBoundaries = length(boundaries);
            bField = bTensor4thOrderField(nBoundaries);
            
            for b = 1:nBoundaries
                
                boundary = boundaries(b);
                kind = boundary.kind;
                radical = 'PatchTensor4thField';
                pField = eval([kind, radical, '(mesh, boundary)']);
                bField = bField.setPatch(b,pField);
            end
        end
    end
end
