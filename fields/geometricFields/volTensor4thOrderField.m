classdef volTensor4thOrderField < tensor4thOrderField
    %volTensor4thOrderField Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % FVM Space-time mesh handle 
        mesh_
        
        % Boundary field
        bTensor4thOrderField_
        
        % Number of interior elements
        nInteriorElements_
        
        % Number of interior + boundary elements
        nTotalElements_
    end
    
    methods
        function obj = volTensor4thOrderField(mesh, input1, input2)
            %volTensor4thOrderField Create a geometric field.
            %
            %   volTensor4thOrderField(mesh) is a geomtric field
            %   with default boundary (calculated/patchTensor4thOrderField).
            %
            %   volTensor4thOrderField(mesh, a_tensor4thOrderField) is the
            %   same as above, except that the field is inialized with
            %   values from a_tensor4thOrderField.
            %
            %   volTensor4thOrderField(mesh, boundaries) is geometric field
            %   with boundary information given by the boundaries array.
            %
            %   volTensor4thOrderField...
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
                
                args{1} = mesh.numberOfElements;
                fromMeshConstructor = true;
            
            % Copy constructor from raw data.
            elseif nargin == 2
                
                % Copy constructor from raw data.
                if isa(input1, 'tensor4thOrderField')
                
                    nElements = mesh.numberOfElements;
                    internalTensorField = ...
                        tensor4thOrderField(input1(1:nElements));
                    args{1} = internalTensorField;
                    copyFromTensorFieldConstructor = true;
                    
                % The second argument must be the boundary conditions set.
                elseif isfield(input1, 'index') ...
                        && isfield(input1, 'kind') ...
                        && length(input1) == mesh.numberOfBoundaries
                    
                    boundaries = input1;
                    args{1} = mesh.numberOfElements;
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
                    volTensor4thOrderField.createDefaultBoundary(obj.mesh_);
                
            elseif copyFromTensorFieldConstructor
                
                obj.bTensor4thOrderField_ = ...
                    volTensor4thOrderField.extractFromTensorField...
                    (...
                        obj.mesh_,...
                        input1...
                    );
            
            elseif fromMeshWithBoundaryInfoConstructor
                
                obj.bTensor4thOrderField_ = ...
                    volTensor4thOrderField.extractBoundaryFieds...
                    (...
                        obj.mesh_,...
                        boundaries...
                    );
                
            elseif fromPiecesConstructor
                    
                obj.bTensor4thOrderField_ = input2;
            else
                error('No possible construction.')
            end
            
            obj.nInteriorElements_ = obj.mesh_.numberOfElements;
            obj.nTotalElements_ = ...
                obj.nInteriorElements_ + obj.mesh_.numberOfBElements;
        end
        
        %% BEG - Data member access
        function mesh = mesh(obj)
            
            mesh = obj.mesh_;
        end
        
        function bTensor4thOrderField = bTensor4thOrderField(obj)
            
            bTensor4thOrderField = obj.bTensor4thOrderField_;
        end
        
        function boundary = boundaryField(obj,b,i,j)
            
            if nargin == 1
                % Hack to support a.boundaryField(1:end)
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
        
        function nInteriorElements = nInteriorElements(obj)
           
            nInteriorElements = obj.nInteriorElements_;
        end
        
        function nTotalElements = nTotalElements(obj)
           
            nTotalElements = obj.nTotalElements_;
        end
        %% END
        
        function output = internalField(obj,i,j)
            
            if nargin == 1

                output = tensor4thOrderField(obj);
                
            elseif nargin == 2

                output = internalField@tensor4thOrderField(obj,i);

            elseif nargin == 3

                output = internalField@tensor4thOrderField(obj,i,j);

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
        
        function newVolTensor4thOrderField = plus(lhs, rhs)
            
            if isa(lhs, 'volTensor4thOrderField') ...
                    && isa(rhs, 'volTensor4thOrderField')
                
                newInternalField = lhs.plus@tensor4thOrderField(rhs);
                
                newBoundaryField = lhs.bTensor4thOrderField_...
                    + rhs.bTensor4thOrderField_;
                
                newVolTensor4thOrderField = ...
                    volTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                errorOperation('Cannot sum operands.');
            end
        end
        
        function newVolTensor4thOrderField = minus(lhs, rhs)
            
            if isa(lhs, 'volTensor4thOrderField') ...
                    && isa(rhs, 'volTensor4thOrderField')
                
                newInternalField = lhs.minus@tensor4thOrderField(rhs);
                
                 newBoundaryField = lhs.bTensor4thOrderField_ ...
                     - rhs.bTensor4thOrderField_;
                
                newVolTensor4thOrderField = ...
                    volTensor4thOrderField...
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
                    volTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newInternalField = rhs.mtimes@tensor4thOrderField(lhs);
                
                newBoundaryField = rhs.bTensor4thOrderField_ * lhs;
                
                newGeometricField = ...
                    volTensor4thOrderField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'volScalarField')

                newInternalField = ...
                    rhs.mtimes@tensor4thOrderField(lhs.internalField());
                
                newBoundaryField = rhs.bTensor4thOrderField_ ...
                    * lhs.boundaryField();
                
                newGeometricField = ...
                    volTensor4thOrderField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(rhs, 'volScalarField')

                newInternalField = ...
                    lhs.mtimes@tensor4thOrderField(rhs.internalField());
                
                newBoundaryField = lhs.bTensor4thOrderField_ ...
                    * rhs.boundaryField();
                
                newGeometricField = ...
                    volTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                errorOperation('Cannot multiply operands.');
            end
        end
        
        function newVolTensor4thOrderField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.mrdivide@tensor4thOrderField(rhs);
                
                newBoundaryField = lhs.bTensor4thOrderField_ / rhs;
                
                newVolTensor4thOrderField = ...
                    volTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'volScalarField')

                newInternalField = ...
                    rhs.mrdivide@tensor4thOrderField(lhs.internalFIeld());
                
                newBoundaryField = rhs.bTensor4thOrderField_ ...
                    / lhs.boundaryField();
                
                newVolTensor4thOrderField = ...
                    volTensor4thOrderField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                errorOperation('Cannot divide operands.');
            end
        end
        
        function newVolumeField = doubleDot(lhs, rhs)
            %doubleDot double constraction like lhs:rhs.
            %   doubleDot(lhs, rhs) returns a volTensorField. This
            %   tensor field being the result of the double contraction:
            %       lhs_{abcd}*rhs_{cd}.
            %   lhs: any volTensor4thOrderField.
            %   rhs: any volTensorField.
            
            if isa(rhs, 'volTensorField')

                newInternalField = ...
                    lhs.doubleDot@tensor4thOrderField(rhs.internalField());
                
                newBoundaryField = lhs.bTensor4thOrderField_ ...
                    .doubleDot(rhs.bTensorField);
                
                newVolumeField = ...
                    volTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                errorOperation('Cannot double dot operands.');
            end
        end
    end
    
    methods (Static)
        
        function bfield = extractFromTensorField(mesh, inputTenField)
            
            nBoundaries = mesh.numberOfBoundaries;
            bfield = bTensor4thOrderField(nBoundaries);
            
            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                startIndex = ...
                    boundary.startFace - mesh.numberOfInteriorFaces...
                    + mesh.numberOfElements;
                nFaces = boundary.numberOfBFaces;
                endIndex = startIndex + nFaces - 1;
                segment = ...
                    tensor4thOrderField(inputTenField(startIndex:endIndex));
                
                pField = patchTensor4thField(mesh, boundary.index, segment);
                bfield = bfield.setPatch(b,pField);
            end
        end
        
        function bField = createDefaultBoundary(mesh)
            
            nBoundaries = mesh.numberOfBoundaries;
            bField = bTensor4thOrderField(nBoundaries);

            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                pField = patchTensor4thField(mesh, boundary.index);
                bField = bField.setPatch(b,pField);
            end
        end
        
        function bField = extractBoundaryFieds(mesh, boundaries) %#ok<INUSL>

            nBoundaries = length(boundaries);
            bField = bTensor4thOrderField(nBoundaries);
            
            for b = 1:nBoundaries
                
                boundary = boundaries(b);
                kind = boundary.kind;
                radical = 'PatchTensor4thOrderField';
                pField = eval([kind, radical, '(mesh, boundary)']);
                bField = bField.setPatch(b,pField);
            end
        end
    end
end
