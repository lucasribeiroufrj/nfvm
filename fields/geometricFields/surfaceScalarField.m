classdef surfaceScalarField
    %surfaceScalarField Summary of this class goes here
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
        
        function obj = surfaceScalarField(mesh, input1, input2)
            %surfaceScalarField Create a geometric field.
            %
            %   surfaceScalarField(mesh) is a geomtric field
            %   with default boundary (calculated/patchScalarField).
            %
            %   surfaceScalarField(mesh, a_scalarField) is the same as
            %   above, except that the field is inialized with values from
            %   a_scalarField.
            %
            %   surfaceScalarField(mesh, boundaries) is geometric field with
            %   boundary information given by the boundaries array.
            %
            %   surfaceScalarField(mesh, a_scalarField, a_bScalarField) 
            %   used to "mount" another field using the "pieces" 
            %   a_scalarField and a_bScalarField. This is the constructor
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
                
                if isa(input1, 'scalarField')

                    assert(input1.numberOfElements == mesh.numberOfFaces);
                    args{1} = input1.range(1:mesh.numberOfInteriorFaces());
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
                
            elseif nargin == 3 && isa(input1, 'scalarField') ...
                    && isa(input2, 'bScalarField')
                
                nInteriorFaces = mesh.numberOfInteriorFaces;
                assert(input1.numberOfElements == nInteriorFaces)
                args{1} = input1;
                fromPiecesConstructor = true;
            else
                error('No possible construction.')
            end
            
            obj.internalField_ = scalarField(args{:});
            obj.mesh_ = mesh;
            
            if fromMeshConstructor == true
                
                obj.boundaryField_ = ...
                    surfaceScalarField.createDefaultBoundary(obj.mesh_);
                    %surfaceScalarField.extractBoundaryFieds(obj.mesh_);
                
            elseif copyFromRawConstructor
                
                obj.boundaryField_ = ...
                    surfaceScalarField.extractFromScalarField...
                    (...
                        obj.mesh_,...
                        input1...
                    );
                
            elseif fromMeshWithBoundaryInfoConstructor
                
                obj.boundaryField_ = ...
                    surfaceScalarField.extractBoundaryFields...
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
        
        function newSurfaceScalarField = plus(lhs, rhs)
            
            if isa(lhs, 'surfaceScalarField') ...
                    && isa(rhs, 'surfaceScalarField')
                
                newInternalField = lhs.internalField_ + rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ + rhs.boundaryField_;
                
                newSurfaceScalarField = ...
                    surfaceScalarField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceScalarField:plus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newSurfaceScalarField = minus(lhs, rhs)
            
            if isa(lhs, 'surfaceScalarField') ...
                    && isa(rhs, 'surfaceScalarField')
                
                newInternalField = ...
                    lhs.internalField_.minus(rhs.internalField_);
                
                newBoundaryField = lhs.boundaryField_ - rhs.boundaryField_;
                
                newSurfaceScalarField = ...
                    surfaceScalarField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceScalarField:plus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newSurfaceField = mtimes(lhs, rhs)
            
            if length(lhs) == 1 && isfloat(lhs)
                
                newInternalField = rhs.internalField_.mtimes(lhs);
                
                newBoundaryField = rhs.boundaryField_ * lhs;
                
                newSurfaceField = ...
                    surfaceScalarField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            
            elseif length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.internalField_.mtimes(rhs);
                
                newBoundaryField = lhs.boundaryField_ * rhs;
                
                newSurfaceField = ...
                    surfaceScalarField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(rhs, 'surfaceTensorField')
                
                newSurfaceField = rhs.mtimes(lhs);
            
            elseif isa(rhs, 'surfaceVectorField')
                
                newSurfaceField = rhs.mtimes(lhs);
            else
                ME = MException('surfaceScalarField:mtimes',...
                        'Cannot mutiply elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newSurfaceScalarField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newInternalField = lhs.internalField_.mrdivide(rhs);
                
                newBoundaryField = lhs.boundaryField_ / rhs;
                
                newSurfaceScalarField = ...
                    surfaceScalarField...
                    (...
                        lhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newInternalField = rhs.internalField_.mrdivide(lhs);
                
                newBoundaryField = lhs / rhs.boundaryField_;
                
                newSurfaceScalarField = ...
                    surfaceScalarField...
                    (...
                        rhs.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
                
            elseif isa(lhs, 'surfaceScalarField') ...
                    && isa(rhs, 'surfaceScalarField')
                
                newInternalField = lhs.internalField_ / rhs.internalField_;
                
                newBoundaryField = lhs.boundaryField_ / rhs.boundaryField_;
                
                newSurfaceScalarField = ...
                    surfaceScalarField...
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
        
        function newSurfaceScalarField = mpower(obj, exponent)

            if length(exponent) == 1 && isfloat(exponent)
                
                newInternalField = obj.internalField_.mpower(exponent);
                
                newBoundaryField = obj.boundaryField_.mpower(exponent);
                
                newSurfaceScalarField = ...
                    surfaceScalarField...
                    (...
                        obj.mesh_,...
                        newInternalField,...
                        newBoundaryField...
                    );
            else
                ME = MException('surfaceScalarField:mpower',...
                        'Cannot compute power. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newSurfaceScalarField = inv(obj)
            
            newInternalField = obj.internalField_.inv();

            newBoundaryField = obj.boundaryField_.inv();

            newSurfaceScalarField = ...
                surfaceScalarField...
                (...
                    obj.mesh_,...
                    newInternalField,...
                    newBoundaryField...
                );
        end
    end
    
    methods(Static)
        
        function bfield = extractFromScalarField(mesh, inputScaField)
            
            nBoundaries = mesh.numberOfBoundaries;
            bfield = bScalarField(nBoundaries);
            
            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                startIndex = boundary.startFace;
                nFaces = boundary.numberOfBFaces;
                endIndex = startIndex + nFaces - 1;
                segment = inputScaField.subField(startIndex:endIndex);
                
                pField = patchScalarField(mesh, boundary.index, segment);
                bfield = bfield.setPatch(b,pField);
            end
        end
        
        function bField = createDefaultBoundary(mesh)
            
            nBoundaries = mesh.numberOfBoundaries;
            bField = bScalarField(nBoundaries);

            for b = 1:nBoundaries

                boundary = mesh.boundaries(b);
                pField = patchScalarField(mesh, boundary.index);
                bField = bField.setPatch(b,pField);
            end
        end
        
        function bField = extractBoundaryFieds(mesh, boundaries) %#ok<INUSL>

            nBoundaries = length(boundaries);
            bField = bScalarField(nBoundaries);
            
            for b = 1:nBoundaries
                
                boundary = boundaries(b);
                kind = boundary.kind;
                radical = 'PatchScalarField';
                pField = eval([kind, radical, '(mesh, boundary)']);
                bField = bField.setPatch(b,pField);
            end
        end
        
        function result = ones(size) %#ok<STOUT,INUSD>
            
            error('Implement me!');
        end
    end
end

