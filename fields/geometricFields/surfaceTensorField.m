classdef surfaceTensorField
    %surfaceTensorField Summary of this class goes here
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
        function obj = surfaceTensorField(mesh, input1, input2)
            %surfaceTensorField Create a geometric field.
            %
            %   surfaceTensorField(mesh) is a geomtric field
            %   with default boundary (calculated/patchTensorField).
            %
            %   surfaceTensorField(mesh, a_tensorField) is the same as
            %   above, except that the field is inialized with values from
            %   a_tensorField.
            %
            %   surfaceTensorField(mesh, boundaries) is geometric field with
            %   boundary information given by the boundaries array.
            %
            %   surfaceTensorField(mesh, a_tensorField, a_bTensorField) 
            %   used to "mount" another field using the "pieces" 
            %   a_tensorField and a_bTensorField. This is the constructor
            %   used when defining overload operator.
            %
            
            fromMeshConstructor = false;
            copyFromRawConstructor = false;
            fromMeshWithBoundaryInfoConstructor = false;
            fromPiecesConstructor = false;
            
            % Contruct from mesh.
            if nargin == 1
                
                args{1} = mesh.numberOfInteriorFaces;
                fromMeshConstructor = true;
            
            % Copy constructor from raw data.
            elseif nargin == 2
                
                if isa(input1, 'tensorField')
                
                    assert(input1.numberOfElements == mesh.numberOfFaces);
                    args{1} = input1.range(1:mesh.numberOfInteriorFaces);
                    copyFromRawConstructor = true;
                    
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
                
            elseif nargin == 3 && isa(input1, 'tensorField') ...
                    && isa(input2, 'bTensorField')
                
                nInteriorFaces = mesh.numberOfInteriorFaces;
                assert(input1.numberOfElements == nInteriorFaces)
                args{1} = input1;
                fromPiecesConstructor = true;
            else
                error('No possible construction.')
            end
                        
            obj.internalField_ = tensorField(args{:});
            obj.mesh_ = mesh;
            
            if fromMeshConstructor == true
                
                obj.boundaryField_ = ...
                    surfaceTensorField.createDefaultBoundary(obj.mesh_);
                    %surfaceTensorField.extractBoundaryFieds(obj.mesh_);
                
            elseif copyFromRawConstructor
                
                obj.boundaryField_ = ...
                    surfaceTensorField.extractFromTensorField...
                    (...
                        obj.mesh_,...
                        input1...
                    );
                
            elseif fromMeshWithBoundaryInfoConstructor
                
                obj.boundaryField_ = ...
                    surfaceTensorField.extractBoundaryFieds...
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
        
        function newSurfaceTensorField = plus(lhs, rhs)
            
            if isa(lhs, 'surfaceTensorField') ...
                    && isa(rhs, 'surfaceTensorField')
                
                newInternalField = lhs.internalField_ + rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ + rhs.boundaryField_;
                
                newSurfaceTensorField = ...
                    surfaceTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceTensorField:plus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newSurfaceTensorField = minus(lhs, rhs)
            
            if isa(lhs, 'surfaceTensorField') ...
                    && isa(rhs, 'surfaceTensorField')
                
                newInternalField = lhs.internalField_ - rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ - rhs.boundaryField_;
                
                newSurfaceTensorField = ...
                    surfaceTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceTensorField:minus',...
                        'Cannot subtract elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newGeometricField = mtimes(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.internalField_ * rhs;
                
                newBoundaryField = lhs.boundaryField_ * rhs;
                
                newGeometricField = ...
                    surfaceTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newInternalField = rhs.internalField_ * lhs;
                
                newBoundaryField = rhs.boundaryField_ * lhs;
                
                newGeometricField = ...
                    surfaceTensorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(rhs, 'surfaceTensorField') ...
                    && isa(lhs, 'surfaceTensorField')
                    
                newInternalField = lhs.internalField_ * rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ * rhs.boundaryField_;
                
                newGeometricField = ...
                    surfaceTensorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'surfaceTensorField') ...
                    && isa(rhs, 'surfaceVectorField')
                
                newInternalField = lhs.internalField_ * rhs.internalField();
                
                newBoundaryField = lhs.boundaryField_ * rhs.boundaryField();
                
                newGeometricField = ...
                    surfaceVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'surfaceScalarField') ...
                    && isa(rhs, 'surfaceTensorField')

                newInternalField = rhs.internalField_ * lhs.internalField();
                
                newBoundaryField = rhs.boundaryField_ * lhs.boundaryField();
                
                newGeometricField = ...
                    surfaceTensorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'surfaceTensorField') ...
                    && isa(rhs, 'surfaceScalarField')

                newInternalField = lhs.internalField_ * rhs.internalField();
                
                newBoundaryField = lhs.boundaryField_ * rhs.boundaryField();
                
                newGeometricField = ...
                    surfaceTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceTensorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newSurfaceTensorField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.internalField_ / rhs;
                
                newBoundaryField = lhs.boundaryField_ / rhs;
                
                newSurfaceTensorField = ...
                    surfaceTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'volScalarField')

                newInternalField = rhs.internalField_ / lhs.internalField();
                
                newBoundaryField = rhs.boundaryField_ / lhs.boundaryField();
                
                newSurfaceTensorField = ...
                    surfaceTensorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceTensorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newSurfaceTensorField = transpose(obj)
            
            newInternalField = obj.internalField_.transpose();
            
            newBoundaryField = obj.boundaryField_.transpose();
                
            newSurfaceTensorField = ...
                surfaceTensorField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
    
        function newSurfaceScalarField = trace(obj)
            
            newInternalField = obj.internalField_.trace();
                
            newBoundaryField = obj.boundaryField_.trace();

            newSurfaceScalarField = ...
                surfaceScalarField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newSurfaceTensorField = inv(obj)
            
            newInternalField = obj.internalField_.inv();
                
            newBoundaryField = obj.boundaryField_.inv();

            newSurfaceTensorField = ...
                surfaceTensorField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newSurfaceScalarField = det(obj)
            
            newInternalField = obj.internalField_.det();
                
            newBoundaryField = obj.boundaryField_.det();

            newSurfaceScalarField = ...
                surfaceScalarField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newSurfaceVectorField = rowAsVector(obj, idx)
            %rowAsVector Return a row as a surfaceVectorField.
            %   rowAsVector(idx) returns the idx-th row as a
            %   surfaceVectorField.
            
            newInternalField = obj.internalField_.rowAsVector(idx);
            
            newBoundaryField = obj.boundaryField_.rowAsVector(idx);

            newSurfaceVectorField = ...
                surfaceVectorField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newSurfaceVectorField = colAsVector(obj, idx)
            %colAsVector Return a column as a surfaceVectorField.
            %   colAsVector(idx) returns the idx-th column as a
            %   surfaceVectorField.
            
            newInternalField = obj.internalField_.colAsVector(idx);
            
            newBoundaryField = obj.boundaryField_.colAsVector(idx);
            
            newSurfaceVectorField = ...
                surfaceVectorField...
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
                surfaceTensorField...
                (...
                    obj.mesh_,...
                    eigVectorsInternalField,...
                    eigVectorsBoundaryField...
                );
            
            eigValues = ...
                surfaceTensorField...
                (...
                    obj.mesh_,...
                    eigValuesInternalField,...
                    eigValuesBoundaryField...
                );
        end
        
        function Identity = identity(obj)
            
            Identity = obj.mesh_.surfaceI();
        end
        
        function obj = copyBoundaryFrom(obj, other)
            %copyBoundaryFrom Assumes the boundaies are compatible.
            
            obj.boundaryField_ = other.boundaryField_;
        end
    end
    
    methods(Static)
       
        function bfield = extractFromTensorField(mesh, inputTensorField)
            
            nBoundaries = mesh.numberOfBoundaries;
            bfield = bTensorField(nBoundaries);
            
            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                startIndex = boundary.startFace;
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
        
        function bField = extractBoundaryFields(mesh, boundaries) %#ok<INUSL>

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
