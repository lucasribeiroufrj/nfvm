classdef scalarField < fields.AbstractScalarField
    %scalarField It is essentially an ordered lists of scalars.
    %   Detailed explanation goes here
    
    properties
        % 1xN array of numbers
        data
    end
    
    methods
        
        function obj = scalarField(input)
            % scalarField An ordered lists of scalars.
            %   scalarField() Creates an uninitilised field. No memory
            %   allocation is done.
            %
            %   scalarField(N) Creates a field with N elements
            %
            %   scalarField(f) Creates a copy of f if f is a scalarField.
            %
            %   scalarField(Array) Creates a field with the content of
            %   Array ( = 1xN array of numbers  ). Basically it converts
            %   MATLAB's array to scalarField. 
            %   Ex.:
            %       Array = zeros(1,10);
            %       ... fill Array ...
            %       s = scalarField(Array);
            
            % Not an empty field
            if nargin ~= 0
              
                % Create a field of size 'input'.
                if length(input) == 1 ...
                   && (isfloat(input) || isinteger(input))
                    
                    obj.data = zeros(1, input);

                % Copy constructor.
                elseif isa(input, 'scalarField')
                    
                    obj.data = input.data;

                % Copy from raw data.
                elseif size(input, 1) == 1 && length(size(input)) == 2
                    
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
        
        function value = mtimes(lhs, rhs)
            
            if isa(rhs, 'tensorField') ...
                    || isa(rhs, 'vectorField')
                
                value = rhs.mtimes(lhs);
                
            elseif isa(lhs, 'scalarField') && isa(rhs, 'scalarField')
                
                value = scalarField(lhs.data .* rhs.data);
                
            elseif length(rhs) == 1 && isfloat(rhs)
                
                value = scalarField(lhs.data * rhs);
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                value = scalarField(lhs * rhs.data);
            else
                ME = MException('scalarField:mtimes',...
                        'Cannot mutiply elements, not implemented.');
                throw(ME);
            end
        end
        
        function value = plus(lhs, rhs)
            
            if isa(lhs, 'fields.AbstractScalarField') ...
                    && isa(rhs, 'fields.AbstractScalarField')
                
                value = scalarField(lhs.data + rhs.data);
            else
                ME = MException('scalarField:plus',...
                        ['Cannot sum elements.', ...
                        ' Are the inputs fields.AbstractScalarField?']);
                throw(ME);
            end
        end
        
        function value = minus(lhs, rhs)
            
            if isa(lhs, 'scalarField') && isa(rhs, 'scalarField')
                
                value = scalarField(lhs.data - rhs.data);
            else
                ME = MException('scalarField:mtimes',...
                        'Cannot subtract elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function value = mrdivide(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                value = scalarField(lhs.data / rhs);
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                value = scalarField(lhs ./ rhs.data);
                
            elseif isa(lhs, 'scalarField') && isa(rhs, 'scalarField')
                
                value = scalarField(lhs.data ./ rhs.data);
            else
                ME = MException('scalarField:mrdivide',...
                        'Cannot divide elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function value = uminus(obj)
           
            value = scalarField(-obj.data);
        end
        
        function value = mpower(obj, exponent)
            
            newData = obj.data.^exponent;

            value = scalarField(newData);
        end
        
        function value = inv(obj)
           
            newData = 1./obj.data;
            
            value = scalarField(newData);
        end
        
        function value = sqrt_(obj)
            
            value = scalarField.sqrt(obj);
        end
        
        function value = log_(obj)
            %log Natural logarithm.
            %   log() Returns a scalar field of the field.
           
            value = scalarField.log(obj);
        end
        
        function endVal = end(obj, k, n)
            
            if n == 1
                endVal = size(obj.data, 2);
            else
                ME = MException('scalarField:end',...
                        'Not implemented yet.');
                throw(ME);
            end
        end
        
        function nElements = numberOfElements(obj)
            
            nElements = size(obj.data, 2);
        end
        
        function output = internalField(obj, i)

            if nargin == 2

                output = obj.data(i);
            else
                error('Use only one parameter.');
            end
        end
        
        function output = range(obj, range)

            if nargin == 2

                output = obj.data(range);
            else
                error('Use one parameter.');
            end
        end
        
        function output = subField(obj, range)

            output = scalarField(obj.range(range));
        end
        
        function output = field(obj, i)

            if nargin == 2

                output = obj.data(i);
            else
                error('Use only one parameter.');
            end
        end
        
        function ones = ones(obj)
           
            ones = scalarField.ones_(obj.numberOfElements());
        end
        
        function array = double(anyScalarField)
           
            array = anyScalarField.data;
        end
    end
    
    methods(Static)
        
        function value = ones_(nElements)
            
            persistent hashMap;
            
            if isempty(hashMap)
               hashMap = ...
                   containers.Map('KeyType', 'single', 'ValueType', 'any');
            end
            
            if isKey(hashMap, nElements)
                value = hashMap(nElements);
            else
                value = scalarField(ones(1, nElements));
                hashMap(nElements) = value;
            end
        end
        
        function value = sqrt(field)
           
            value = scalarField(sqrt(field.data));
        end
        
        function value = log(field)
           
            value = scalarField(log(field.data));
        end
        
        function value = test(field)
           
            value = scalarField(log(field.data));
        end
    end
    
end

