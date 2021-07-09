classdef bScalarField
    %bScalarField Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % One boundary field for each patch
        patchFields_
        
        nPatches_;
    end
    
    methods
        function obj = bScalarField(nPatches)
            %bScalarField Construct an instance of this class
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
                
            elseif ~isa(field,'patchScalarField')
                
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
            
            newBoundaryField = bScalarField(nPatches);
            
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
            
            newBoundaryField = bScalarField(nPatches);
            
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
                
                newBoundaryField = bScalarField(nPatches);

                for b = 1:nPatches

                    boundary = lhs.patchFields_{b};

                    newBoundary =  boundary * rhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(rhs, 'bVectorField')
                
                newBoundaryField = rhs.mtimes(lhs);
                
            else
                ME = MException('volScalarField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newBoundaryField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bScalarField(nPatches);

                for b = 1:nPatches

                    boundary = lhs.patchFields_{b};

                    newBoundary =  boundary / rhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                nPatches = rhs.nPatches_;
                
                newBoundaryField = bScalarField(nPatches);

                for b = 1:nPatches

                    boundary = rhs.patchFields_{b};

                    newBoundary =  lhs / boundary;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(lhs,'bScalarField')...
                    && isa(rhs,'bScalarField')
                
                nPatches = length(lhs.patchFields_);
            
                if nPatches ~= length(rhs.patchFields_)
                    error('Number of boundaries does not match');
                end

                newBoundaryField = bScalarField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};
                    boundary2 = rhs.patchFields_{b};

                    newBoundary =  boundary1 / boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
            else
                ME = MException('volScalarField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newBoundaryField = mpower(obj, exponent)
            
            if length(exponent) == 1 && isfloat(exponent)
                
                nPatches = obj.nPatches_;
                
                newBoundaryField = bScalarField(nPatches);

                for b = 1:nPatches

                    boundary = obj.patchFields_{b};

                    newBoundary =  boundary.mpower(exponent);
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
            else
                ME = MException('bScalarField:mpower',...
                        'Cannot compute power. Not implemented yet.');
                throw(ME);
            end
        end
        
        function newBoundaryField = inv(obj)
            
            nPatches = obj.nPatches_;

            newBoundaryField = bScalarField(nPatches);

            for b = 1:nPatches

                boundary = obj.patchFields_{b};

                newBoundary =  boundary.inv();
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function value = norm(obj) %#ok<STOUT,MANU>
           
            error('Implement me!');
        end
        
        function obj = updateDivImpKgradU(obj, solidModel)
           
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
        
        function nTotalElements = nTotalElements(obj)
           
            nTotalElements = 0;
            nPatches = obj.nPatches_;
            
            for b = 1:nPatches
                
                nTotalElements = nTotalElements + ...
                    obj.patchFields_{b}.numberOfFaces();
            end
        end
    end
end

