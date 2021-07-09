classdef bTensor4thOrderField
    %bTensor4thOrderField Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % One boundary field for each patch
        patchTensor4thOrderFields_
        
        nPatches_;
    end
    
    methods
        function obj = bTensor4thOrderField(nPatches)
            %bTensor4thOrderField Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.nPatches_ = nPatches;
        end
        
        %% BEG - Data member access
        function obj = setPatch(obj, index, field)
            
            if index > obj.nPatches_
                
                error('Out-of-bound access');
                
            elseif ~isa(field,'patchTensor4thOrderField')
                
                error('Only patch fields can be set.');
            end
            
            obj.patchTensor4thOrderFields_{index} = field;
        end
        
        function nPatches = nPatches(obj)
            
            nPatches = obj.nPatches_;
        end
        %% END
        
        function boundary = boundary(obj, index)
            
            if index > obj.nPatches_
                error('Out-of-bound access');
            end
            
            boundary = obj.patchTensor4thOrderFields_{index};
        end
        
        function newBoundaryField = plus(lhs, rhs)
            
            nPatches = length(lhs.patchTensor4thOrderFields_);
            
            if nPatches ~= length(rhs.patchTensor4thOrderFields_)
                error('Number of boundaries does not match');
            end
            
            newBoundaryField = bTensor4thOrderField(nPatches);
            
            for b = 1:nPatches
                
                boundary1 = lhs.patchTensor4thOrderFields_{b};
                boundary2 = rhs.patchTensor4thOrderFields_{b};

                newBoundary =  boundary1 + boundary2;
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = minus(lhs, rhs)
            
            nPatches = length(lhs.patchTensor4thOrderFields_);
            
            if nPatches ~= length(rhs.patchTensor4thOrderFields_)
                error('Number of boundaries does not match');
            end
            
            newBoundaryField = bTensor4thOrderField(nPatches);
            
            for b = 1:nPatches
                
                boundary1 = lhs.patchTensor4thOrderFields_{b};
                boundary2 = rhs.patchTensor4thOrderFields_{b};

                newBoundary =  boundary1 - boundary2;
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = mtimes(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bTensor4thOrderField(nPatches);

                for b = 1:nPatches

                    boundary = lhs.patchTensor4thOrderFields_{b};

                    newBoundary =  boundary * rhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                nPatches = rhs.nPatches_;
                
                newBoundaryField = bTensor4thOrderField(nPatches);

                for b = 1:nPatches

                    boundary = rhs.patchTensor4thOrderFields_{b};

                    newBoundary =  boundary * lhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(lhs, 'bScalarField')
                    
                nPatches = rhs.nPatches_;
                
                newBoundaryField = bTensor4thOrderField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchTensor4thOrderFields_{b};
                    boundary2 = rhs.patchFields(b);

                    newBoundary =  boundary1 * boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(rhs, 'bScalarField')

                nPatches = lhs.nPatches_;
                
                newBoundaryField = bTensor4thOrderField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchTensor4thOrderFields_{b};
                    boundary2 = rhs.patchFields(b);

                    newBoundary =  boundary1 * boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
            else
                errorOperation('Cannot multiply operands.');
            end
        end
        
        function newBoundaryField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bTensor4thOrderField(nPatches);

                for b = 1:nPatches

                    boundary = lhs.patchTensor4thOrderFields_{b};

                    newBoundary =  boundary / rhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(rhs, 'bScalarField')
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bTensor4thOrderField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchTensor4thOrderFields_{b};
                    boundary2 = lhs.patchFields(b);

                    newBoundary =  boundary1 / boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
            else
                errorOperation('Cannot divide operands.');
            end
        end
        
        function newBoundaryField = doubleDot(lhs, rhs)
            
            if isa(rhs, 'bTensorField')
                
                nPatches = lhs.nPatches_;

                newBoundaryField = bTensorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchTensor4thOrderFields_{b};
                    boundary2 = rhs.patchFields(b);

                    newBoundary =  boundary1.doubleDot(boundary2);
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
            else
                errorOperation('Cannot double dot operands.');
            end
        end
        
        function obj = correctBoundaryConditions(obj, volU, volGradU)
            
            for b = 1:obj.nPatches_
                
                obj = ...
                    obj.setPatch...
                    (...
                        b,...
                        obj.patchTensor4thOrderFields_{b}.correct...
                        (...
                            volU,...
                            volGradU...
                        )...
                    );
            end
        end
    end
end

