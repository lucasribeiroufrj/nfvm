classdef patchTensorField
    %patchTensorField An abstraction of a boundary condition.
    %   Detailed explanation goes here
    %
    % TODO:
    % - We should use only a sub-mesh, but should also have a ref. to the
    % main mesh.
    
    properties
        
        tensorField_;
        userName_;
        index_;
        numberOfBFaces_;
        startFace_;
        endFace_;
        mesh_;
    end
    
    methods
        
        function obj = patchTensorField(mesh, index, inputTensorField)
            %patchTensorField Construct an instance of this class
            %
            %   patchTensorField(mesh, index) is a patch with geometric
            %   boundary information given by
            %   mesh.bondaries(index).
            %
            %   patchTensorField(mesh, index, inputVectorField) is the same 
            %   as above with addition that the field is initialized with
            %   inputTensorField.
            %
            
            boundary = mesh.boundaries(index);
            nFaces = boundary.numberOfBFaces;
            
            % Copy from a boundary.
            if nargin == 2
               
                args{1} = nFaces;
                
            % Initialize the boundary with values from the tensorField.
            elseif nargin == 3 && isa(inputTensorField, 'tensorField')
                
                assert(inputTensorField.numberOfElements == nFaces);
                args{1} = inputTensorField;
            else
                error('No possible construction.')
            end
            
            obj.tensorField_ = tensorField(args{:});
            
            obj.userName_ = boundary.userName;
            obj.index_ = boundary.index;
            obj.numberOfBFaces_ = nFaces;
            obj.startFace_ = boundary.startFace;
            obj.endFace_ = obj.startFace_ + obj.numberOfBFaces_ - 1;
            obj.mesh_ = mesh;
        end
        
        %% BEG - Data member access
        function field = tensorField(obj)
            
            field = obj.tensorField_;
        end
        
        function obj = setTensorField(obj, newValue)
            
            obj.tensorField_ = newValue;
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
            
                output = obj.tensorField_;
                
            elseif nargin == 2

                output = obj.tensorField_.field(i);

            elseif nargin == 3

                output = obj.tensorField_.field(i,j);

            else
                error('Use one or two parameters');
            end
        end
        
        function newPatch = plus(lhs, rhs)
            
            if isa(lhs, 'patchTensorField') ...
                    && isa(rhs, 'patchTensorField')
                
                newTensorField = ...
                    lhs.tensorField_.plus(rhs.tensorField_);
                
                newPatch = ...
                    patchTensorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = minus(lhs, rhs)
           
            if isa(lhs, 'patchTensorField') ...
                    && isa(rhs, 'patchTensorField')
                
                newTensorField = ...
                    lhs.tensorField_.minus(rhs.tensorField_);
                
                newPatch = ...
                    patchTensorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = mtimes(lhs, rhs)
           
            if length(rhs) == 1 && isfloat(rhs)
                
                newTensorField = lhs.tensorField_.mtimes(rhs);
                
                newPatch = ...
                    patchTensorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newTensorField = rhs.tensorField_.mtimes(lhs);
                
                newPatch = ...
                    patchTensorField...
                    (...
                        rhs.mesh_,...
                        rhs.index_,...
                        newTensorField...
                    );
                
            elseif isa(lhs, 'patchTensorField')
                
                if isa(rhs, 'patchTensorField')
                
                    newTensorField = ...
                        lhs.tensorField_.mtimes(rhs.tensorField_);

                    newPatch = ...
                        patchTensorField...
                        (...
                            lhs.mesh_,...
                            lhs.index_,...
                            newTensorField...
                        );
                
                elseif isa(rhs, 'patchScalarField')
                
                    newTensorField = ...
                        lhs.tensorField_.mtimes(rhs.field());

                    newPatch = ...
                        patchTensorField...
                        (...
                            lhs.mesh_,...
                            lhs.index_,...
                            newTensorField...
                        );
                elseif isa(rhs, 'patchVectorField')
                
                    newVectorField = ...
                        lhs.tensorField_.mtimes(rhs.field());

                    newPatch = ...
                        patchVectorField...
                        (...
                            rhs.mesh_,...
                            rhs.index_,...
                            newVectorField...
                        );
                end
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = mrdivide(lhs, rhs)
           
            if length(rhs) == 1 && isfloat(rhs)
                
                newTensorField = lhs.tensorField_.mrdivide(rhs);
                
                newPatch = ...
                    patchTensorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
                
            elseif isa(rhs, 'patchScalarField')
                
                newTensorField = ...
                    lhs.tensorField_.mrdivide(rhs.field());
                
                newPatch = ...
                    patchTensorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
            else
                error('Not implemented yet');
            end
        end
        
        function newPatch = transpose(obj)
            
            newTensorField = obj.tensorField_.transpose();
                
            newPatch = ...
                patchTensorField...
                (...
                    obj.mesh_,...
                    obj.index_,...
                    newTensorField...
                );
        end
        
        function newPatch = trace(obj)
            
            newScalarField = obj.tensorField_.trace();
                
            newPatch = ...
                patchScalarField...
                (...
                    obj.mesh_,...
                    obj.index_,...
                    newScalarField...
                );
        end
        
        function newPatch = inv(obj)
            
            newTensorField = obj.tensorField_.inv();
                
            newPatch = ...
                patchTensorField...
                (...
                    obj.mesh_,...
                    obj.index_,...
                    newTensorField...
                );
        end
        
        function newPatch = det(obj)
            
            newScalarField = obj.tensorField_.det();
                
            newPatch = ...
                patchScalarField...
                (...
                    obj.mesh_,...
                    obj.index_,...
                    newScalarField...
                );
        end
        
        function newPatch = rowAsVector(obj, idx)
            %rowAsVector Return a row as a patchVectorField.
            %   rowAsVector(idx) returns the idx-th row as a
            %   patchVectorField.
            
            newVectorField = obj.tensorField_.rowAsVector(idx);
                
            newPatch = ...
                patchVectorField...
                (...
                    obj.mesh_,...
                    obj.index_,...
                    newVectorField...
                );
        end
        
        function newPatch = colAsVector(obj, idx)
            %colAsVector Return a column as a patchVectorField.
            %   colAsVector(idx) returns the idx-th column as a
            %   patchVectorField.
            
            newVectorField = obj.tensorField_.colAsVector(idx);
                
            newPatch = ...
                patchVectorField...
                (...
                    obj.mesh_,...
                    obj.index_,...
                    newVectorField...
                );
        end
        
        function [pEigVectorsField, pEigValuesField] = eigen(obj)
            
            [ newEigVectorsTensorField, newEigValuesTensorField] = ...
                obj.tensorField_.eigen();
                
            pEigVectorsField = ...
                patchTensorField...
                (...
                    obj.mesh_,...
                    obj.index_,...
                    newEigVectorsTensorField...
                );
            
            pEigValuesField = ...
                patchTensorField...
                (...
                    obj.mesh_,...
                    obj.index_,...
                    newEigValuesTensorField...
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
        
        function Identity = identity(obj)
            
            Identity = obj.mesh_.volI().boundaryField(obj.index());
        end
    end
end

