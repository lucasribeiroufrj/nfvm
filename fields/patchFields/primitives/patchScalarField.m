classdef patchScalarField
    %patchScalarField An abstraction of a boundary condition.
    %   Detailed explanation goes here
    %
    % TODO:
    % - We should use only a sub-mesh, but should also have a ref. to the
    % main mesh.
    
    properties
        
        scalarField_;
        userName_;
        index_;
        numberOfBFaces_;
        startFace_;
        endFace_;
        mesh_;
    end
    
    methods
        
        function obj = patchScalarField(mesh, index, inputScalarField)
            %patchScalarField Construct an instance of this class
            %
            %   patchScalarField(mesh, index) is a patch with geometric
            %   boundary information given by
            %   mesh.bondaries(index).
            %
            %   patchScalarField(mesh, index,inputScalarField) is the same 
            %   as above with addition that the field is initialized with
            %   inputScalarField.
            %
            
            boundary = mesh.boundaries(index);
            nFaces = boundary.numberOfBFaces;
            
            % Copy from a boundary.
            if nargin == 2

                args{1} = nFaces;
                
            % Initialize the boundary with values from the scalarField.
            elseif nargin == 3 && isa(inputScalarField, 'scalarField')
                
                assert(inputScalarField.numberOfElements == nFaces);
                args{1} = inputScalarField;
            else
                error('No possible construction.')
            end
            
            obj.scalarField_ = scalarField(args{:});
            
            obj.userName_ = boundary.userName;
            obj.index_ = boundary.index;
            obj.numberOfBFaces_ = nFaces;
            obj.startFace_ = boundary.startFace;
            obj.endFace_ = obj.startFace_ + obj.numberOfBFaces_ - 1;
            obj.mesh_ = mesh;
        end
        
        %% BEG - Data member access
        function field = scalarField(obj)
            
            field = obj.scalarField_;
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
            
                output = obj.scalarField_;
                
            elseif nargin == 2

                output = obj.scalarField_.field(i);

            elseif nargin == 3

                output = obj.scalarField_.field(i,j);

            else
                error('Use one or two parameters');
            end
        end
        
        function newPatch = plus(lhs, rhs)
            
            if isa(lhs, 'patchScalarField') ...
                    && isa(rhs, 'patchScalarField')
                
                newScalarField = ...
                    lhs.scalarField_.plus(rhs.scalarField_);
                
                newPatch = ...
                    patchScalarField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newScalarField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = minus(lhs, rhs)
           
            if isa(lhs, 'patchScalarField') ...
                    && isa(rhs, 'patchScalarField')
                
                newScalarField = ...
                    lhs.scalarField_.minus(rhs.scalarField_);
                
                newPatch = ...
                    patchScalarField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newScalarField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = mtimes(lhs, rhs)
           
            if length(rhs) == 1 && isfloat(rhs)
                
                newScalarField = lhs.scalarField_.mtimes(rhs);
                
                newPatch = ...
                    patchScalarField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newScalarField...
                    );
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newScalarField = rhs.scalarField_.mtimes(lhs);
                
                newPatch = ...
                    patchScalarField...
                    (...
                        rhs.mesh_,...
                        rhs.index_,...
                        newScalarField...
                    );
                
            elseif isa(rhs, 'patchTensorField') ||  ...
                isa(rhs, 'patchVectorField')
                
                newPatch = rhs.mtimes(lhs);
                
            elseif isa(lhs, 'patchScalarField') ...
                    && isa(rhs, 'patchScalarField')

            newScalarField = ...
                lhs.scalarField_.mtimes(rhs.field());

            newPatch = ...
                patchScalarField...
                (...
                    lhs.mesh_,...
                    lhs.index_,...
                    newScalarField...
                );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                newScalarField = lhs.scalarField_.mrdivide(rhs);
                
                newPatch = ...
                    patchScalarField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newScalarField...
                    );
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newScalarField = rhs.scalarField_.mrdivide(lhs);
                
                newPatch = ...
                    patchScalarField...
                    (...
                        rhs.mesh_,...
                        rhs.index_,...
                        newScalarField...
                    );
                
            elseif isa(lhs, 'patchScalarField') ...
                    && isa(rhs, 'patchScalarField')
                
                newScalarField = ...
                    lhs.scalarField_.mrdivide(rhs.field());
                
                newPatch = ...
                    patchScalarField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newScalarField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = mpower(obj, exponent)
           
            if length(exponent) == 1 && isfloat(exponent)
                
                newScalarField = obj.scalarField_.mpower(exponent);
                
                newPatch = ...
                    patchScalarField...
                    (...
                        obj.mesh_,...
                        obj.index_,...
                        newScalarField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = inv(obj)
                
            newScalarField = obj.scalarField_.inv();

            newPatch = ...
                patchScalarField...
                (...
                    obj.mesh_,...
                    obj.index_,...
                    newScalarField...
                );
        end
        
        function value = norm(obj) %#ok<STOUT,MANU>
           
            error('Implement me!');
        end
        
        function nTotalElements = nTotalElements(obj)
           
            nTotalElements = obj.numberOfBFaces_;
        end
        
        function obj = updateDivImpKgradU(obj, solidModel) %#ok<INUSD>
            % updateDivImpKgradU Computes the div(K*gradU).
            
            error('not implemented');
        end
        
        function obj = correct(obj, volU, volGradU)
            
            error('not implemented');
        end
    end
    
    methods(Static)
        
        function result = ones(size) %#ok<STOUT,INUSD>
            
            error('Implement me!');
        end
    end
end

