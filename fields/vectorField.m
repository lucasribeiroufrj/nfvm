classdef vectorField
    %vectorField It is essentially an ordered lists of vectors.
    %   Detailed explanation goes here
    
    properties
        % 3xN array of numbers
        data
    end
    
    methods
        
        function obj = vectorField(input)
            % vectorField An ordered lists of vectors.
            %   vectorField() Creates an uninitilised field. No memory
            %   allocation is done.
            %
            %   vectorField(N) Creates a field with N elements
            %
            %   vectorField(f) Creates a copy of f if f is a vectorField.
            %
            %   vectorField(Array) Creates a field with the content of
            %   Array ( = 3xN array of numbers  ). Basically it converts
            %   MATLAB's array to vectorField. 
            %   Ex.:
            %       Array = zeros(3,10);
            %       ... fill Array ...
            %       s = vectorField(Array);
            
            % Not an empty field
            if nargin ~= 0
                
                % Create a field of size 'input'.
                if length(input) == 1 ...
                   && (isfloat(input) || isinteger(input))
                    
                    obj.data = zeros(3, input);

                % Copy constructor.
                elseif isa(input, 'vectorField')
                    
                    obj.data = input.data;

                % Copy from raw data.
                elseif size(input, 1) == 3 && length(size(input)) == 2
                    
                    obj.data = input;
                else
                    ME = MException('vectorField:constructor',...
                        'Invalid parameter');
                    throw(ME);
                end
            else
                % This is an empty field, no need to initialize it.
            end
        end
        
        function data = getData(obj)
           
            data = obj.data;
        end
        
        function obj = setData(obj, newValue)
           
            obj.data = newValue;
        end
        
        function value = vector_X_scalarPar(lhs,rhs)
            
            newData = ...
                vectorField.mulVecxScaLinear...
                (...
                    lhs.data(1,:),...
                    lhs.data(2,:),...
                    lhs.data(3,:),...
                    rhs.data(:) ...
                );
            
            value = vectorField(newData);
        end
        
        function value = vector_and_vector(lhs,rhs)
           
            newData = ...
                reshape...
                (...
                    vectorField.andVecxVecLinear...
                    (...
                        lhs.data(1,:),...
                        lhs.data(2,:),...
                        lhs.data(3,:),...
                        rhs.data(1,:),...
                        rhs.data(2,:),...
                        rhs.data(3,:) ...
                    ),3,3,lhs.numberOfElements()...
                );
            
            value = tensorField(newData);
        end
        
        function value = mtimes(lhs, rhs)
            
            global VECTORIZE
            
            if length(rhs) == 1 && isfloat(rhs)
                
                value = vectorField(lhs.data * rhs);
            
            elseif length(lhs) == 1 && isfloat(lhs)
                
                value = vectorField(rhs.data * lhs);
            
            elseif isa(rhs, 'scalarField')
                
                if VECTORIZE
                    value = vector_X_scalarPar(lhs,rhs);
                else
                    value = ...
                        vectorField(bsxfun(@mtimes, lhs.data, rhs.data));
                end
            else
                ME = MException('vectorField:mtimes',...
                        'Cannot mutiply elements');
                throw(ME);
            end
        end
        
        function value = and(lhs, rhs)
            %and Tensorial product between two vectors.
            
            global VECTORIZE
            
            if isa(lhs, 'vectorField') && isa(lhs, 'vectorField') 
                
                if VECTORIZE
                    value = vector_and_vector(lhs,rhs);
                else
                    N = size(lhs.data,2);
                    rawData = zeros(3,3,N);
                    ldata = lhs.data;
                    rdata = rhs.data;

                    for i=1:N

                        rawData(:,:,i) = ldata(:,i) * rdata(:,i).';
                    end

                    value = tensorField(rawData);
                end
            else
                ME = MException('vectorField:and',...
                        'Cannot mutiply elements');
                throw(ME);
            end
        end
        
        function value = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                value = vectorField(lhs.data / rhs);
            
            elseif length(lhs) == 1 && isfloat(lhs)
                error('Impossible');
                %value = vectorField(rhs.data / lhs);
                
            elseif isa(rhs, 'scalarField')
                
                value = vectorField(bsxfun(@rdivide,lhs.data,rhs.data));
            else
                ME = MException('vectorField:mrdivide',...
                        'Cannot divide elements');
                throw(ME);
            end
        end
        
        function value = plus(lhs, rhs)
            
            if isa(lhs, 'vectorField') && isa(rhs, 'vectorField')
                
                value = vectorField(lhs.data + rhs.data);
            else
                ME = MException('vectorField:mtimes',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function value = minus(lhs, rhs)
            
            if isa(lhs, 'vectorField') && isa(rhs, 'vectorField')
                
                value = vectorField(lhs.data - rhs.data);
            else
                ME = MException('vectorField:mtimes',...
                        'Cannot subtract elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function value = uminus(obj)
           
            value = vectorField(-obj.data);
        end
        
        function value = mpower(lhs, rhs)
            
            if isa(lhs, 'vectorField') && isa(rhs, 'vectorField')
                
                value = scalarField(dot(lhs.data, rhs.data));
            else
                ME = MException('vectorField:mpower',...
                        'Cannot calculate the inner product.');
                throw(ME);
            end
        end
        
        function value = norm(obj)
            
            p = 2; % Euclidean norm
            direction = 1; % Column norm
            value = sqrt(sum(abs(obj.data).^p, direction));
            value = scalarField(value);
        end
        
        function endVal = end(obj, k, n)
            
            if n == 1
                endVal = size(obj.data, 2);
            else
                ME = MException('vectorField:end',...
                        'Not implemented yet.');
                throw(ME);
                %endVal = builtin('end', obj.data, k, n);
            end
        end
        
        function nElements = numberOfElements(obj)
            
            nElements = size(obj.data, 2);
        end
        
        function output = internalField(obj, i, j)

            if nargin == 2

                output = obj.data(:,i);

            elseif nargin == 3

                output = obj.data(i,j);

            else
                error('Use one or two parameters');
            end
        end
        
        function output = range(obj, range)

            if nargin == 2

                output = obj.data(:,range);
            else
                error('Use one parameter.');
            end
        end
        
        function output = subField(obj, range)

            output = vectorField(obj.range(range));
        end
        
        function output = field(obj, i, j)

            if nargin == 2

                output = obj.data(:,i);

            elseif nargin == 3

                output = obj.data(i,j);

            else
                error('Use one or two parameters');
            end
        end
        
        function obj = setFromRawData(obj, data, fstIdx, sndIdx)
            % setFromRawData(data) return the same field with the new data.
            %   data must accept the indices (1:3, fstIdx:sndIdx) and
            %   fstIdx will be mapped onto the element of index 1.
            
            if fstIdx > sndIdx || fstIdx < 1
                error('invalid indexing');
            end
            
            nElements = obj.numberOfElements();
            
            if (sndIdx - fstIdx + 1) > nElements
                
                error('invalid indexing');
            end
            
            nElements = sndIdx - fstIdx + 1;
            obj.data(1:3, 1:nElements) = data(1:3, fstIdx:sndIdx);
        end
        
        function i = i(obj)
            %i returns the vectorField x such that x(i) = [1 0 0].
           
            i = vectorField.i_(obj.numberOfElements());
        end
        
        function j = j(obj)
            %j returns the vectorField x such that x(i) = [0 1 0].
           
            j = vectorField.j_(obj.numberOfElements());
        end
        
        function k = k(obj)
            %k returns the vectorField x such that x(i) = [0 0 1].
           
            k = vectorField.k_(obj.numberOfElements());
        end
    end
    
    methods (Static)
        
        function i = i_(nElements)
            %i_ returns a vectorField x such that x(i) = [1 0 0]
            %   i_(nElements) nElements must be the size of the desired
            %   field x.
            
            i = vectorField.base(1, nElements);
        end
        
        function j = j_(nElements)
            %j_ returns a vectorField x such that x(i) = [0 1 0]
            %   i_(nElements) nElements must be the size of the desired
            %   field x.
            
            j = vectorField.base(2, nElements);
        end
        
        function k = k_(nElements)
            %k_ returns a vectorField x such that x(i) = [0 0 1]
            %   i_(nElements) nElements must be the size of the desired
            %   field x.
            
            k = vectorField.base(3, nElements);
        end
       
        function x = base(dir,nElements)
            %base returns a vectorField x such that x(i) = [1 0 0] if
            %   dir = 1, x(i) = [0, 1, 0] if dir = 2 and etc.
            %   base(nElements) nElements must be the size of the desired
            %   field x.
            
            persistent hashMaps;
            
            if isempty(hashMaps)
               hashMaps = cell(3,1);
            end
            
            if isempty(hashMaps{dir})
                hashMaps{dir} = ...
                    containers.Map('KeyType','single','ValueType','any');
            end
            
            if isKey(hashMaps{dir}, nElements)
                x = hashMaps{dir}(nElements);
            else
                data = zeros([3,nElements]);
                data(dir,:) = 1;
                x = vectorField(data);
                hashMaps{dir}(nElements) = x;
            end
        end
        
        function output = mulVecxScaLinear(v1,v2,v3,a)
            
            output = ...
                [...
                    v1.*a';...
                    v2.*a';...
                    v3.*a'; ...
                ];
        end
        
        function output = andVecxVecLinear(v1,v2,v3,u1,u2,u3)
            
            output = ...
                [...
                    u1.*v1;u1.*v2;u1.*v3;...
                    u2.*v1;u2.*v2;u2.*v3;...
                    u3.*v1;u3.*v2;u3.*v3 ...
                ];
        end
    end
end

