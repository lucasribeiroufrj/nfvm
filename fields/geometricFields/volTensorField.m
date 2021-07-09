classdef volTensorField
    %volTensorField Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Interior field
        internalField_
        
        % Boundary field
        boundaryField_
        
        % FVM Space-time mesh handle 
        mesh_
    end
    
    methods
        function obj = volTensorField(mesh, input1, input2)
            %volTensorField Create a geometric field.
            %
            %   volTensorField(mesh) is a geomtric field
            %   with default boundary (calculated/patchScalarField).
            %
            %   volTensorField(mesh, a_tensorField) is the same as
            %   above, except that the field is inialized with values from
            %   a_tensorField.
            %
            %   volTensorField(mesh, boundaries) is geometric field with
            %   boundary information given by the boundaries array.
            %
            %   volTensorField(mesh, a_tensorField, a_bTensorField) used to
            %   "mount" another field using the "pieces" a_tensorField and
            %   a_bTensorField. This is the constructor used when defining
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
            
            % Copy constructor from raw data.
            elseif nargin == 2
                
                % Copy constructor from raw data.
                if isa(input1, 'tensorField')
                
                    nTotalElements = mesh.numberOfTotalElements();
                    assert(input1.numberOfElements == nTotalElements);
                    args{1} = input1.range(1:mesh.numberOfElements());
                    copyFromRawConstructor = true;
                    
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
                
            elseif nargin == 3 && isa(input1, 'tensorField') ...
                    && isa(input2, 'bTensorField')
                
                nElements = mesh.numberOfElements();
                assert(input1.numberOfElements() == nElements);
                args{1} = input1;
                fromPiecesConstructor = true;
            else
                error('No possible construction.')
            end
            
            obj.internalField_ = tensorField(args{:});
            obj.mesh_ = mesh;
            
            if fromMeshConstructor == true
                
                obj.boundaryField_ = ...
                    volTensorField.createDefaultBoundary(obj.mesh_);
                    %volTensorField.extractBoundaryFieds(obj.mesh_);
                
            elseif copyFromRawConstructor
                
                obj.boundaryField_ = ...
                    volTensorField.extractFromTensorField...
                    (...
                        obj.mesh_,...
                        input1...
                    );
            
            elseif fromMeshWithBoundaryInfoConstructor
                
                obj.boundaryField_ = ...
                    volTensorField.extractBoundaryFieds...
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
        
        function bTensorField = bTensorField(obj)
            
            bTensorField = obj.boundaryField_;
        end
        
        function output = internalField(obj,i,j)
            
            if nargin == 1
                % Hack to support a.internalField(1:end)
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
            % field Internal usage function
            %   It is to support, e.g., a.boundaryField(2).field(1:2)
            %
            
            value = obj;
        end
        
        function newVolTensorField = plus(lhs, rhs)
            
            if isa(lhs, 'volTensorField') && isa(rhs, 'volTensorField')
                
                newInternalField = lhs.internalField_ + rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ + rhs.boundaryField_;
                
                newVolTensorField = ...
                    volTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volTensorField:plus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newVolTensorField = minus(lhs, rhs)
            
            if isa(lhs, 'volTensorField') && isa(rhs, 'volTensorField')
                
                newInternalField = lhs.internalField_ - rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ - rhs.boundaryField_;
                
                newVolTensorField = ...
                    volTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volTensorField:minus',...
                        'Cannot subtract elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newGeometricField = mtimes(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.internalField_ * rhs;
                
                newBoundaryField = lhs.boundaryField_ * rhs;
                
                newGeometricField = ...
                    volTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newInternalField = rhs.internalField_ * lhs;
                
                newBoundaryField = rhs.boundaryField_ * lhs;
                
                newGeometricField = ...
                    volTensorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'volTensorField') && isa(rhs, 'volTensorField')
                    
                newInternalField = lhs.internalField_ * rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ * rhs.boundaryField_;
                
                newGeometricField = ...
                    volTensorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'volTensorField') && isa(rhs, 'volVectorField')
                
                newInternalField = lhs.internalField_ * rhs.internalField();
                
                newBoundaryField = lhs.boundaryField_ * rhs.boundaryField();
                
                newGeometricField = ...
                    volVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'volScalarField') && isa(rhs, 'volTensorField')

                newInternalField = rhs.internalField_ * lhs.internalField();
                
                newBoundaryField = rhs.boundaryField_ * lhs.boundaryField();
                
                newGeometricField = ...
                    volTensorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'volTensorField') && isa(rhs, 'volScalarField')

                newInternalField = lhs.boundaryField_ * rhs.internalField();
                
                newBoundaryField = lhs.boundaryField_ * rhs.boundaryField();
                
                newGeometricField = ...
                    volTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volTensorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newVolTensorField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.internalField_ / rhs;
                
                newBoundaryField = lhs.boundaryField_ / rhs;
                
                newVolTensorField = ...
                    volTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'volScalarField')

                newInternalField = rhs.internalField_ / lhs.internalField();
                
                newBoundaryField = rhs.boundaryField_ / lhs.boundaryField();
                
                newVolTensorField = ...
                    volTensorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('volTensorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newVolTensorField = transpose(obj)
            
            newInternalField = obj.internalField_.transpose();
            
            newBoundaryField = obj.boundaryField_.transpose();
                
            newVolTensorField = ...
                volTensorField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newVolScalarField = trace(obj)
            
            newInternalField = obj.internalField_.trace();
                
            newBoundaryField = obj.boundaryField_.trace();

            newVolScalarField = ...
                volScalarField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newVolTensorField = inv(obj)
            
            newInternalField = obj.internalField_.inv();
                
            newBoundaryField = obj.boundaryField_.inv();

            newVolTensorField = ...
                volTensorField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newVolScalarField = det(obj)
            
            newInternalField = obj.internalField_.det();
                
            newBoundaryField = obj.boundaryField_.det();

            newVolScalarField = ...
                volScalarField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newVolVectorField = rowAsVector(obj, idx)
            %rowAsVector Return a row as a volVectorField.
            %   rowAsVector(idx) returns the idx-th row as a
            %   volVectorField.
            
            newInternalField = obj.internalField_.rowAsVector(idx);
            
            newBoundaryField = obj.boundaryField_.rowAsVector(idx);

            newVolVectorField = ...
                volVectorField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newVolVectorField = colAsVector(obj, idx)
            %colAsVector Return a column as a volVectorField.
            %   colAsVector(idx) returns the idx-th column as a
            %   volVectorField.
            
            newInternalField = obj.internalField_.colAsVector(idx);
            
            newBoundaryField = obj.boundaryField_.colAsVector(idx);
            
            newVolVectorField = ...
                volVectorField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function [eigVectors, eigValues] = eigen(obj)
            
            [eigVectorsInternalField, eigValuesInternalField] = ...
                obj.internalField_.eigen();

            [eigVectorsBoundaryField , eigValuesBoundaryField] = ...
                obj.boundaryField_.eigen();

            eigVectors = ...
                volTensorField...
                (...
                    obj.mesh_,...
                    eigVectorsInternalField,...
                    eigVectorsBoundaryField...
                );
            
            eigValues = ...
                volTensorField...
                (...
                    obj.mesh_,...
                    eigValuesInternalField,...
                    eigValuesBoundaryField...
                );
        end
        
        function Identity = identity(obj)
            
            Identity = obj.mesh_.volI();
        end
        
        function obj = copyBoundaryFrom(obj, other)
            %copyBoundaryFrom Assumes the boundaies are compatible.
            
            obj.boundaryField_ = other.boundaryField_;
        end
    end
    
    methods (Static)
        
        function bfield = extractFromTensorField(mesh, inputTensorField)
            
            nBoundaries = mesh.numberOfBoundaries;
            bfield = bTensorField(nBoundaries);
            
            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                startIndex = ...
                    boundary.startFace - mesh.numberOfInteriorFaces...
                    + mesh.numberOfElements;
                nFaces = boundary.numberOfBFaces;
                endIndex = startIndex + nFaces - 1;
                segment = inputTensorField.subField(startIndex:endIndex);
                
                pField = patchTensorField(mesh, boundary.index, segment);
                bfield = bfield.setPatch(b,pField);
            end
        end
        
        function bField = createDefaultBoundary(mesh)
            
            nBoundaries = mesh.numberOfBoundaries;
            bField = bTensorField(nBoundaries);

            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                pField = patchTensorField(mesh, boundary.index);
                bField = bField.setPatch(b,pField);
            end
        end
        
        function bField = extractBoundaryFieds(mesh, boundaries) %#ok<INUSL>

            nBoundaries = length(boundaries);
            bField = bTensorField(nBoundaries);
            
            for b = 1:nBoundaries
                
                boundary = boundaries(b);
                kind = boundary.kind;
                radical = 'PatchTensorField';
                pField = eval([kind, radical, '(mesh, boundary)']);
                bField = bField.setPatch(b,pField);
            end
        end
    end
end
