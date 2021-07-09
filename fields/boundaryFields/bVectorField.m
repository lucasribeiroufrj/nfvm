classdef bVectorField
    %bVectorField Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % One boundary field for each patch
        patchFields_;
        
        nPatches_;
    end
    
    methods
        function obj = bVectorField(nPatches)
            %bVectorField Construct an instance of this class
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
                
            elseif ~isa(field,'patchVectorField')
                
                error('Only patch fields can be setPatched.');
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
            
            newBoundaryField = bVectorField(nPatches);
            
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
            
            newBoundaryField = bVectorField(nPatches);
            
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
                
                newBoundaryField = bVectorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};

                    newBoundary =  boundary1 * rhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(rhs,'bScalarField')
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bVectorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};
                    boundary2 = rhs.patchField(b);

                    newBoundary =  boundary1 * boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            else
                ME = MException('volVectorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newBoundaryField = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bVectorField(nPatches);

                for b = 1:nPatches

                    boundary = lhs.patchFields_{b};

                    newBoundary =  boundary / rhs;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            elseif isa(rhs, 'bScalarField')
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bVectorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};
                    boundary2 = rhs.patchField(b);

                    newBoundary =  boundary1 / boundary2;

                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
                
            else
                ME = MException('volVectorField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function newBoundaryField = mpower(lhs, rhs)
            
            if isa(lhs, 'bVectorField') && isa(rhs, 'bVectorField')
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bScalarField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};
                    boundary2 = rhs.patchFields_{b};

                    newBoundary =  boundary1 ^ boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
            else
                ME = MException('volVectorField:mpower',...
                        'Cannot calculate the inner product.');
                throw(ME);
            end
        end
        
        function newBoundaryField = and(lhs, rhs)
            
            if isa(lhs, 'bVectorField') && isa(rhs, 'bVectorField')
                
                nPatches = lhs.nPatches_;
                
                newBoundaryField = bTensorField(nPatches);

                for b = 1:nPatches

                    boundary1 = lhs.patchFields_{b};
                    boundary2 = rhs.patchFields_{b};

                    newBoundary =  boundary1 & boundary2;
                    newBoundaryField = ...
                        newBoundaryField.setPatch(b,newBoundary);
                end
            else
                ME = MException('volVectorField:and',...
                        'Cannot calculate the dyadic product.');
                throw(ME);
            end
        end
        
        function newBoundaryField = uminus(obj)
           
            nPatches = obj.nPatches_;
                
            newBoundaryField = bVectorField(nPatches);

            for b = 1:nPatches

                boundary = obj.patchFields_{b};
                newBoundary =  boundary.uminus();
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function newBoundaryField = norm(obj)
           
            nPatches = obj.nPatches_;
                
            newBoundaryField = bScalarField(nPatches);

            for b = 1:nPatches

                boundary = obj.patchFields_{b};
                newBoundary =  boundary.norm();
                newBoundaryField = ...
                    newBoundaryField.setPatch(b,newBoundary);
            end
        end
        
        function obj = updateDivImpKgradUfluxes(obj, solidModel)
            % updateDivImpKgradUfluxes Computes the coeff. of div(K*gradU)
            % term.
            
            for b = 1:obj.nPatches_
                
                updatedPatch = ...
                    obj.patchFields_{b}.updateDivImpKgradUfluxes...
                    (...
                        solidModel...
                    );
                
                obj = ...
                    obj.setPatch...
                    (...
                        b,...
                        updatedPatch...
                    );
            end
        end
        
        function obj = ...
                correctBoundaryConditions...
                (...
                    obj,...
                    volU,...
                    volGradU,...
                    surfaceGradU...
                )
            
            for b = 1:obj.nPatches_
                
                obj = ...
                    obj.setPatch...
                    (...
                        b,...
                        obj.patchFields_{b}.correct...
                        (...
                            volU,...
                            volGradU,...
                            surfaceGradU...
                        )...
                    );
            end
        end
        
        function obj = setFromRawData(obj, data, fstIdx, sndIdx)
            % setFromRawData(data) return the same field with the new data.
            %   data must accept the indices (1:3, fstIdx:sndIdx) and
            %   fstIdx will be mapped onto the boundary field face of index
            %   1.
            
            if fstIdx > sndIdx || fstIdx < 1
                error('invalid indexing');
            end
            
            nElements = 0;
            nPatches = obj.nPatches_;
            
            for b = 1:nPatches
                
                nElements = nElements + ...
                    obj.patchFields_{b}.numberOfFaces();
            end
            
            if (sndIdx - fstIdx + 1) > nElements
                error('invalid indexing');
            end

            % Each patch "consumes" its own parcel.
            for b = 1:nPatches

                nBFaces = obj.patchFields_{b}.numberOfFaces();
                rgtIdx = fstIdx + nBFaces - 1;
                
                if fstIdx <= rgtIdx
                
                    if sndIdx <= rgtIdx

                        k = sndIdx;
                        left = 0;
                    else
                        left = sndIdx - rgtIdx;
                        k = rgtIdx;
                    end

                    obj.patchFields_{b} = obj. ...
                        patchFields_{b}.setFromRawData...
                        (...
                            data,...
                            fstIdx,...
                            k...
                        );
                    
                    fstIdx = k + 1;
                else
                    left = sndIdx - rgtIdx;
                end
                
                if left <= 0
                    return;
                end
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

