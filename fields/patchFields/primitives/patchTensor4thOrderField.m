classdef patchTensor4thOrderField < tensor4thOrderField
    %patchTensor4thOrderField An abstraction of a boundary condition.
    %   Detailed explanation goes here
    %
    % TODO:
    % - We should use only a sub-mesh, but should also have a ref. to the
    % main mesh.
    
    properties
        
        userName_;
        index_;
        numberOfBFaces_;
        startFace_;
        endFace_;
        mesh_;
    end
    
    methods
        
        function obj = ...
                patchTensor4thOrderField(mesh, index, inputTensorField)
            %patchTensor4thOrderField Construct an instance of this class
            %
            %   patchTensor4thOrderField(mesh, index) is a patch with
            %   geometric boundary information given by
            %   mesh.bondaries(index).
            %
            %   patchTensor4thOrderField...
            %   (...
            %       mesh,...
            %       index,...
            %       inputTensor4thOrderField...
            %   )
            %   is the same as above with addition that the field is
            %   initialized with inputTensor4thOrderField.
            %
            
            boundary = mesh.boundaries(index);
            nFaces = boundary.numberOfBFaces;
            
            % Copy from a boundary.
            if nargin == 2
               
                args{1} = nFaces;
                
            % Initialize the boundary with values from the
            % tensor4thOrderField.
            elseif nargin == 3 ...
                    && isa(inputTensorField, 'tensor4thOrderField')
                
                args{1} = inputTensorField(1:nFaces);
            else
                error('No possible construction.')
            end
            
            obj@tensor4thOrderField(args{:});
            
            obj.userName_ = boundary.userName;
            obj.index_ = boundary.index;
            obj.numberOfBFaces_ = nFaces;
            obj.startFace_ = boundary.startFace;
            obj.endFace_ = obj.startFace_ + obj.numberOfBFaces_ - 1;
            obj.mesh_ = mesh;
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
        
        function output = field(obj,i,j)
            
            if nargin == 1
            
                output = tensor4thOrderField(obj);
                
            elseif nargin == 2

                output = field@tensor4thOrderField(obj,i);

            elseif nargin == 3

                output = field@tensor4thOrderField(obj,i,j);

            else
                error('Use one or two parameters');
            end
        end
        
        function newPatch = plus(lhs, rhs)
            
            if isa(lhs, 'patchTensor4thOrderField') ...
                    && isa(rhs, 'patchTensor4thOrderField')
                
                newTensorField = lhs.plus@tensor4thOrderField(rhs);
                
                newPatch = ...
                    patchTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
            else
                errorOperation('Cannot subtract operands.');
            end
        end
        
        function newPatch = minus(lhs, rhs)
           
            if isa(lhs, 'patchTensor4thOrderField')...
                    && isa(rhs, 'patchTensor4thOrderField')
                
                newTensorField = lhs.minus@tensor4thOrderField(rhs);
                
                newPatch = ...
                    patchTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
            else
                errorOperation('Cannot subtract operands.');
            end
        end
        
        function newPatch = mtimes(lhs, rhs)
           
            if length(rhs) == 1 && isfloat(rhs)
                
                newTensorField = lhs.mtimes@tensor4thOrderField(rhs);
                
                newPatch = ...
                    patchTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                newTensorField = rhs.mtimes@tensor4thOrderField(lhs);
                
                newPatch = ...
                    patchTensor4thOrderField...
                    (...
                        rhs.mesh_,...
                        rhs.index_,...
                        newTensorField...
                    );
                
            elseif isa(lhs, 'patchTensor4thOrderField') ...
                    && isa(rhs, 'patchScalarField')
                
                newTensorField = lhs.mtimes@tensor4thOrderField(rhs);
                
                newPatch = ...
                    patchTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
                
            elseif isa(lhs, 'patchVectorField') ...
                    && isa(rhs, 'patchTensor4thOrderField')
                
                newVectorField = rhs.mtimes@tensor4thOrderField(lhs);
                
                newPatch = ...
                    patchVectorField...
                    (...
                        rhs.mesh_,...
                        rhs.index_,...
                        newVectorField...
                    );
            else
                errorOperation('Cannot multiply operands.');
            end
        end
        
        function newPatch = mrdivide(lhs, rhs)
           
            if length(rhs) == 1 && isfloat(rhs)
                
                newTensorField = lhs.mrdivide@tensor4thOrderField(rhs);
                
                newPatch = ...
                    patchTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
                
            elseif isa(rhs, 'patchScalarField')
                
                newTensorField = lhs.mrdivide@tensor4thOrderField(rhs);
                
                newPatch = ...
                    patchTensor4thOrderField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
            else
                errorOperation('Cannot divide operands.');
            end
        end
        
        function newPatch = doubleDot(lhs, rhs)
            %doubleDot double constraction like lhs:rhs.
            %   doubleDot(lhs, rhs) returns a patchTensor4thOrderField.
            %   This tensor field being the result of the double
            %   contraction:
            %       lhs_{abcd}*rhs_{cd}.
            %   lhs: any patchTensor4thOrderField.
            %   rhs: any patchTensorField.
           
            if isa(rhs, 'patchTensorField')
                
                newTensorField = lhs.doubleDot@tensor4thOrderField(rhs);
                
                newPatch = ...
                    patchTensorField...
                    (...
                        lhs.mesh_,...
                        lhs.index_,...
                        newTensorField...
                    );
            else
                errorOperation('Cannot double dot operands.');
            end
        end
        
        function mesh = mesh(obj)
            
            mesh = obj.mesh_;
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
end

