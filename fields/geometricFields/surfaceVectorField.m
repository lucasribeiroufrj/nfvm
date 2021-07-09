classdef surfaceVectorField
    %surfaceVectorField Summary of this class goes here
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
        
        function obj = surfaceVectorField(mesh, input1, input2)
            %surfaceVectorField Create a geometric field.
            %
            %   surfaceVectorField(mesh) is a geomtric field
            %   with default boundary (calculated/patchVectorField).
            %
            %   surfaceVectorField(mesh, a_vectorField) is the same as
            %   above, except that the field is inialized with values from
            %   a_vectorField.
            %
            %   surfaceVectorField(mesh, boundaries) is geometric field with
            %   boundary information given by the boundaries array.
            %
            %   surfaceVectorField(mesh, a_vectorField, a_bVectorField) 
            %   used to "mount" another field using the "pieces" 
            %   a_vectorField and a_bVectorField. This is the constructor
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
                
                if isa(input1, 'vectorField')
                
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
                
            elseif nargin == 3 && isa(input1, 'vectorField') ...
                    && isa(input2, 'bVectorField')
                
                nInteriorFaces = mesh.numberOfInteriorFaces;
                assert(input1.numberOfElements == nInteriorFaces)
                args{1} = input1;
                fromPiecesConstructor = true;
            else
                error('No possible construction.')
            end
            
            obj.internalField_ = vectorField(args{:});
            obj.mesh_ = mesh;
            
            if fromMeshConstructor == true
                
                obj.boundaryField_ = ...
                    surfaceVectorField.createDefaultBoundary(obj.mesh_);
                
            elseif copyFromRawConstructor
                
                obj.boundaryField_ = ...
                    surfaceVectorField.extractFromVectorField...
                    (...
                        obj.mesh_,...
                        input1...
                    );
                
            elseif fromMeshWithBoundaryInfoConstructor
                
                obj.boundaryField_ = ...
                    surfaceVectorField.extractBoundaryFieds...
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
        
        function value = field(obj)
            % field internal usage function
            %   It is to support, e.g., a.boundaryField(2).field(1:2)
            %
            
            value = obj;
        end
        
        function newSurfaceVectorField = plus(lhs, rhs)
            
            if isa(lhs, 'surfaceVectorField') ...
                    && isa(rhs, 'surfaceVectorField')
                
                newInternalField = ...
                    lhs.internalField_.plus(rhs.internalField_);
                
                newBoundaryField = ...
                    lhs.boundaryField_ + rhs.boundaryField_;
                
                newSurfaceVectorField = ...
                    surfaceVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceVectorField:plus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newSurfaceVectorField = minus(lhs, rhs)
            
            if isa(lhs, 'surfaceVectorField') && isa(rhs, 'surfaceVectorField')
                
                newInternalField = ...
                    lhs.internalField_.minus(rhs.internalField_);
                
                newBoundaryField = ...
                    lhs.boundaryField_ - rhs.boundaryField_;
                
                newSurfaceVectorField = ...
                    surfaceVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceVectorField:minus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newSurfaceVectorField = mtimes(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.internalField_.mtimes(rhs);
                
                newBoundaryField = lhs.boundaryField_ * rhs;
                
                newSurfaceVectorField = ...
                    surfaceVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif length(lhs) == 1 && isfloat(lhs)

                newInternalField = rhs.internalField_.mtimes(lhs);
                
                newBoundaryField = rhs.boundaryField_ * (lhs);
                
                newSurfaceVectorField = ...
                    surfaceVectorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(rhs, 'surfaceScalarField')
                
                newInternalField = lhs.internalField_ * rhs.internalField();
                
                newBoundaryField = lhs.boundaryField_ * rhs.boundaryField();
                
                newSurfaceVectorField = ...
                    surfaceVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceVectorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newSurfaceVectorField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.internalField_.mrdivide(rhs);
                
                newBoundaryField = lhs.boundaryField_ / rhs;
                
                newSurfaceVectorField = ...
                    surfaceVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newInternalField = rhs.internalField_.mrdivide(lhs);
                
                newBoundaryField = lhs / rhs.boundaryField_;
                
                newSurfaceVectorField = ...
                    surfaceVectorField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(rhs, 'surfaceScalarField')
                
                newInternalField = lhs.internalField_ / rhs.internalField();
                
                newBoundaryField = lhs.boundaryField_ / rhs.boundaryField();
                
                newSurfaceVectorField = ...
                    surfaceVectorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceVectorField:mrdivide',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newSurfaceScalarField = mpower(lhs, rhs)
            
            if isa(lhs, 'surfaceVectorField') ...
                    && isa(rhs, 'surfaceVectorField')
                
                newInternalField = ...
                    lhs.internalField_.mpower(rhs.internalField_);
                
                newBoundaryField = lhs.boundaryField_ ^ rhs.boundaryField_;
                
                newSurfaceScalarField = ...
                    surfaceScalarField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceVectorField:mpower',...
                        'Cannot calculate the inner product.');
                throw(ME);
            end
        end
        
        function newSurfaceTensorField = and(lhs, rhs)
            %and Tensorial product between two vectors.
            
            if isa(lhs, 'surfaceVectorField') ...
                    && isa(rhs, 'surfaceVectorField')
                
                newInternalField = lhs.internalField_ & rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ & rhs.boundaryField_;
                
                newSurfaceTensorField = ...
                    surfaceTensorField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceVectorField:and',...
                        'Cannot calculate the dyadic product.');
                throw(ME);
            end
        end
        
        function newVectorField = uminus(obj)
           
            newInternalField = obj.internalField_.uminus();
                
            newBoundaryField = obj.boundaryField_.uminus();

            newVectorField = ...
                surfaceVectorField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function newSurfaceScalarField = norm(obj)
           
            newInternalField = obj.internalField_.norm();
                
            newBoundaryField = obj.boundaryField_.norm();

            newSurfaceScalarField = ...
                surfaceScalarField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
        
        function obj = updateDivImpKgradU(obj, solidModel)
           
            error('Implement me!');
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
        
        function obj = setInternalData(obj, newValue)
            
            obj.internalField_ = obj.internalField_.setData(newValue);
        end
        
        function obj = setInternalField(obj, newValue)
            
            obj.internalField_ = newValue;
        end
    end
    
    methods (Static)
        
        function bfield = extractFromVectorField(mesh, inputVecField)
            
            nBoundaries = mesh.numberOfBoundaries;
            bfield = bVectorField(nBoundaries);
            
            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                startIndex = boundary.startFace;
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
                
                boundary = boundaries(b);
                kind = boundary.kind;
                radical = 'PatchVectorField';
                pField = eval([kind, radical, '(mesh, boundary)']);
                bField = bField.setPatch(b,pField);
            end
        end
    end
    
end

