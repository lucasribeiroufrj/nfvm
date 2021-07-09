classdef bTensorField
    %bTensorField Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % One boundary field for each patch
        patchFields_
        
        nPatches_;
    end
    
    methods
        function obj = bTensorField(nPatches)
            %bTensorField Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.nPatches_ = nPatches;
        end
        
        %% BEG - Data member access
        function patchFields = patchFields(obj)
            %patchFields Return the patch cells
            
            patchFields = obj.patchFields_;
        end
        
        function nPatches = nPatches(obj)
            %nPatches Return the number of patches
            
            nPatches = obj.nPatches_;
        end
        
        function obj = setPatch(obj, index, field)
            
            if index > obj.nPatches_
                
                error('Out-of-bound access');
                
            elseif ~isa(field,'patchTensorField')
                
                error('Only patch fields can be set.');
            end
            
            obj.patchFields_{index} = field;
        end
        %% END - Data member access
        
        function patchField = patchField(obj, index)
            %patchField Return the index-th patch as a field.
            
            patchField = obj.patchFields_{index};
        end
        
        function boundary = boundary(obj, index)
            
            if index > obj.nPatches_
                error('Out-of-bound access');
            end
            
            boundary = obj.patchFields_{index};
        end
        
        function newBoundaryField = plus(lhs, rhs)
            
            nPatches = length(lhs.patchFields_);
            
            if nPatches ~= length(rhs.patchFields_)
                error('Number of boundaries does not match');
            end
            
            newBoundaryField = bTensorField(nPatches);
            
            for b = 1:nPatches
                
                boundary1 = lhs.patchFields_{b};
                boundary2 = rhs.patchFields_{b};

                newBoundary =  boundary1 + boundary2;
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = minus(lhs, rhs)
            
            nPatches = length(lhs.patchFields_);
            
            if nPatches ~= length(rhs.patchFields_)
                error('Number of boundaries does not match');
            end
            
            newBoundaryField = bTensorField(nPatches);
            
            for b = 1:nPatches
                
                boundary1 = lhs.patchFields_{b};
                boundary2 = rhs.patchFields_{b};

                newBoundary =  boundary1 - boundary2;
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = mtimes(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bTensorField(nPatches);

                for b = 1:nPatches

                    boundary = lhs.patchFields_{b};

                    newBoundary =  boundary * rhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                nPatches = rhs.nPatches_;
                
                newBoundaryField = bTensorField(nPatches);

                for b = 1:nPatches

                    boundary = rhs.patchFields_{b};

                    newBoundary =  boundary * lhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(lhs, 'bTensorField') && isa(rhs, 'bTensorField')
                
                nPatches = rhs.nPatches_;
                
                newBoundaryField = bTensorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};
                    boundary2 = rhs.patchFields_{b};

                    newBoundary =  boundary1 * boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(lhs, 'bTensorField') && isa(rhs, 'bVectorField')
                
                nPatches = rhs.nPatches_;
                
                newBoundaryField = bVectorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};
                    boundary2 = rhs.patchField(b);

                    newBoundary =  boundary1 * boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(lhs, 'bScalarField') && isa(rhs, 'bTensorField')
                    
                nPatches = rhs.nPatches_;
                
                newBoundaryField = bTensorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchField(b);
                    boundary2 = rhs.patchFields_{b};

                    newBoundary =  boundary1 * boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(lhs, 'bTensorField') && isa(rhs, 'bScalarField')

                nPatches = lhs.nPatches_;
                
                newBoundaryField = bTensorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};
                    boundary2 = rhs.patchField(b);

                    newBoundary =  boundary1 * boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
            else
                ME = MException('bTensorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newBoundaryField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bTensorField(nPatches);

                for b = 1:nPatches

                    boundary = lhs.patchFields_{b};

                    newBoundary =  boundary / rhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(rhs, 'bScalarField')
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bTensorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};
                    boundary2 = lhs.patchField(b);

                    newBoundary =  boundary1 / boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
            else
                ME = MException('bTensorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newBoundaryField = transpose(obj)
            
            nPatches = obj.nPatches_;
                
            newBoundaryField = bTensorField(nPatches);

            for b = 1:nPatches

                boundary = obj.patchFields_{b};

                newBoundary =  boundary.transpose();
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = trace(obj)
            
            nPatches = obj.nPatches_;
                
            newBoundaryField = bScalarField(nPatches);

            for b = 1:nPatches

                boundary = obj.patchFields_{b};

                newBoundary =  boundary.trace();
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = inv(obj)
            
            nPatches = obj.nPatches_;
                
            newBoundaryField = bTensorField(nPatches);

            for b = 1:nPatches

                boundary = obj.patchFields_{b};

                newBoundary =  boundary.inv();
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = det(obj)
            
            nPatches = obj.nPatches_;
                
            newBoundaryField = bScalarField(nPatches);

            for b = 1:nPatches

                boundary = obj.patchFields_{b};

                newBoundary =  boundary.det();
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = rowAsVector(obj, idx)
            %rowAsVector Return a row as a bVectorField.
            %   rowAsVector(idx) returns the idx-th row as a
            %   bVectorField.
            
            nPatches = obj.nPatches_;
                
            newBoundaryField = bVectorField(nPatches);

            for b = 1:nPatches

                boundary = obj.patchFields_{b};

                newBoundary =  boundary.rowAsVector(idx);
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = colAsVector(obj, idx)
            %colAsVector Return a column as a bVectorField.
            %   colAsVector(idx) returns the idx-th column as a
            %   bVectorField.
            
            nPatches = obj.nPatches_;
                
            newBoundaryField = bVectorField(nPatches);

            for b = 1:nPatches

                boundary = obj.patchFields_{b};

                newBoundary =  boundary.colAsVector(idx);
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function [bEigVectorsField, bEigValuesField] = eigen(obj)
            
            nPatches = obj.nPatches_;
                
            bEigVectorsField = bTensorField(nPatches);
            bEigValuesField = bTensorField(nPatches);
            
            for b = 1:nPatches
                
                boundary = obj.patchFields_{b};
                [ new_bEigenVectors, new_bEigenValues ] = boundary.eigen();
                
                bEigVectorsField = ...
                    bEigVectorsField.setPatch(b, new_bEigenVectors);
                bEigValuesField = ...
                    bEigValuesField.setPatch(b, new_bEigenValues);
            end
        end
        
        function value = norm(obj) %#ok<STOUT,MANU>
           
            error('Implement me!');
        end
        
        function obj = updateDivImpKgradU(obj, solidModel)
           
            % Create a new patch and setPatch one by one.
            
            for b = 1:obj.nPatches_
                
                obj = ...
                    obj.setPatch...
                    (...
                        b,...
                        obj.patchFields_{b}.updateDivImpKgradU...
                        (...
                            solidModel...
                        )...
                    );
            end
        end
        
        function obj = correctBoundaryConditions(obj, volU, volGradU)
            
            for b = 1:obj.nPatches_
                
                obj = ...
                    obj.setPatch...
                    (...
                        b,...
                        obj.patchFields_{b}.correct...
                        (...
                            volU,...
                            volGradU...
                        )...
                    );
            end
        end
    end
    
    methods (Static)
        
        function Identity = Identity(mesh)
           
            error('Implement me!');
        end
    end
end

