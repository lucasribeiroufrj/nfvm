classdef tensor4thOrderField
    %tensor4thOrderField It is essentially an ordered lists of tensors.
    %   Detailed explanation goes here
    
    properties
        data
    end
    
    methods
        
        function obj = tensor4thOrderField(input)
            % tensor4thOrderField An ordered lists of 4th order tensors.
            %   tensor4thOrderField() Creates an uninitilised field. No
            %   memory allocation is done.
            %
            %   tensor4thOrderField(N) Creates a field with N elements
            %
            %   tensor4thOrderField(f) Creates a copy of f if f is a
            %   tensor4thOrderField.
            %
            %   tensor4thOrderField(Array) Creates a field with the content
            %   of Array ( = 3x3x3x3xN array of numbers  ). Basically it
            %   converts MATLAB's array to tensor4thOrderField. 
            %   Ex.:
            %       Array = zeros(3,3,3,3,10);
            %       ... fill Array ...
            %       s = tensor4thOrderField(Array);
            
            % Not an empty field
            if nargin ~= 0
              
                % Create a field of size 'input'.
                if length(input) == 1 ...
                   && (isfloat(input) || isinteger(input))
                    
                    obj.data = zeros(3,3,3,3,input);

                % Copy constructor.
                elseif isa(input, 'tensor4thOrderField')
                    
                    obj.data = input.data;

                % Copy from raw data.
                elseif size(input, 1) == 3 ...
                        && size(input, 2) == 3 ...
                        && size(input, 3) == 3 ...
                        && size(input, 4) == 3 ...
                        && length(size(input)) == 5
                    
                    obj.data = input;
                else
                    ME = MException('tensor4thOrderField:constructor',...
                        'Invalid parameters.');
                    throw(ME);
                end
            else
                % This is an empty field, no need to initialize it.
            end
        end
        
        function value = plus(lhs, rhs)
            
            if isa(lhs, 'tensor4thOrderField') ...
                    && isa(rhs, 'tensor4thOrderField')
                
                value = tensor4thOrderField(lhs.data + rhs.data);
            else
                errorOperation('Cannot sum operands.');
            end
        end
        
        function value = mtimes(lhs, rhs)
            
            if length(rhs) == 1 && isfloat(rhs)
                
                value = tensor4thOrderField(lhs.data * rhs);
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                value = tensor4thOrderField(lhs * rhs.data);
                
            elseif isa(rhs, 'scalarField')
                
                nElements = rhs.numberOfElements();
                
                newData = lhs.data;
                ldata = lhs.data;
                rdata = rhs.data;
                
                for i=1:nElements
                    newData(:,:,:,:,i) = ldata(:,:,:,:,i) * rdata(i);
                end
                
                value = tensor4thOrderField(newData);
                
            elseif isa(rhs, 'tensorField')
                
                nElements = rhs.numberOfElements();
                
                % Memory allocation.
                newData = lhs.data;
                
                % For performance.
                ldata = lhs.data;
                rdata = rhs.data;
                
                for i=1:nElements
                    for j=1:3
                        newData(:,:,:,j,i) = ...
                              ldata(:,:,:,1,i)*rdata(1,j,i)...
                            + ldata(:,:,:,2,i)*rdata(2,j,i)...
                            + ldata(:,:,:,3,i)*rdata(3,j,i);
                    end
                end
                
                value = tensor4thOrderField(newData);
            else
                errorOperation('Cannot mutiply operands.');
            end
        end
        
        function nElements = numberOfElements(obj)
            
            nElements = size(obj.data, 5);
        end
        
        function newTensorField = doubleDot(lhs, rhs)
            %doubleDot double constraction like lhs:rhs.
            %   doubleDot(lhs, rhs) returns a tensorField. This tensor
            %   field being the result of the double contraction:
            %       lhs_{abcd}*rhs_{cd}.
            %   lhs: any tensor4thOrderField.
            %   rhs: any tensorField.
            
            if isa(rhs, 'tensorField')
                
                nElements = size(lhs.data, 5);

                % Memory Allocation
                newData(3,3,nElements) = rhs.data(1,1,1);

                ldata = lhs.data;
                rdata = rhs.data;

                for i = 1:nElements
                    for a = 1:3
                        for b = 1:3
                            matrix = reshape(ldata(a,b,:,:,i), [3, 3]);
                            newData(a,b,i) = sum(diag(matrix*rdata(:,:,i).'));
                        end
                    end
                end

                newTensorField = tensorField(newData);
            else
                errorOperation('Cannot double dot operands.');
            end
        end
        
        function value = contractWithVector(lhs, rhs, index)
            %contractWithVector Creates tensor3rdOrderField by means of
            %   contraction.
            %
            %   contractWithVector(lhs, rhs, index) returns
            %   "lhs_{abcd}*rhs{a}" if index == 1,
            %   "lhs_{abcd}*rhs{b}" if index == 2, ...
            %   lhs: a tensor4thOrderTensor
            %   lhs: a vectorField
            %   index: any integer from [1 2 3 4]
            
            % Annotation:
            % A = rand([3 3 3])
            % subsref(A,struct('type', '()', 'subs', {{':' ':' [3]}}))
            
            if isa(rhs, 'vectorField')
                tic
                nElements = rhs.numberOfElements();
                
                % Memory allocation.
                newData = zeros(3,3,3,nElements);
                
                % For performance.
                ldata = lhs.data;
                rdata = rhs.data;
                
                if index == 1
                    
                    for i=1:nElements

                        newData(:,:,:,i) = ...
                              ldata(1,:,:,:,i)*rdata(1,i)...
                            + ldata(2,:,:,:,i)*rdata(2,i)...
                            + ldata(3,:,:,:,i)*rdata(3,i);
                    end
                    
                elseif index == 2
                    
                    for i=1:nElements

                        newData(:,:,:,i) = ...
                              ldata(:,1,:,:,i)*rdata(1,i)...
                            + ldata(:,2,:,:,i)*rdata(2,i)...
                            + ldata(:,3,:,:,i)*rdata(3,i);
                    end
                    
                elseif index == 3
                    
                    for i=1:nElements

                        newData(:,:,:,i) = ...
                              ldata(:,:,1,:,i)*rdata(1,i)...
                            + ldata(:,:,2,:,i)*rdata(2,i)...
                            + ldata(:,:,3,:,i)*rdata(3,i);
                    end
                    
                elseif index == 4
                    
                    for i=1:nElements

                        newData(:,:,:,i) = ...
                              ldata(:,:,:,1,i)*rdata(1,i)...
                            + ldata(:,:,:,2,i)*rdata(2,i)...
                            + ldata(:,:,:,3,i)*rdata(3,i);
                    end
                else
                    error('Invalid index.');
                end
                
                value = tensor3rdOrderField(newData);
            else
                errorOperation('Cannot contract operands.');
            end
        end
    end
    
    methods (Static)
       
        function newTensor4thOrderField = compose(lhs, rhs, n)
            %compose Creates a 4th order tensor from two tensor fields.
            %
            %   compose(lhs,rhs,n) returns a tensor4thOrderField.
            %   lhs: a tensorField.
            %   rhs: a tensorField.
            %   n: any permutation of [1 2 3 4].
            %   So, compose(A,B,[1 3 2 4]) := C{1,2,3,4} = A_{1,3}*B_{2,4}.
            
            if isa(lhs, 'tensorField') && isa(rhs, 'tensorField')
                
                nElements = lhs.numberOfElements();
                data = zeros(3,3,3,3,nElements);
                ldata = lhs.data;
                rdata = rhs.data;
                
                for z = 1:nElements
                    for i = 1:3
                        for j = 1:3
                            for k = 1:3
                                for l = 1:3
                                    
                                    if nargin == 2
                                        data(i,j,k,l,z) = ...
                                            ldata(i, j, z)...
                                            * rdata(k, l, z);
                                    else
                                        a = [i j k l];
                                        data(i,j,k,l,z) = ...
                                            ldata(a(n(1)), a(n(2)), z)...
                                            * rdata(a(n(3)), a(n(4)), z);
                                    end
                                end
                            end
                        end
                    end
                end
                
                newTensor4thOrderField = tensor4thOrderField(data);
            else
                errorOperation('Cannot make the tensor product.');
            end
        end
        
        function newTensorField = doubleDotSym(lhs, rhs)
            %doubleDotSym double constraction like lhs:rhs.
            %   doubleDotSym(lhs, rhs) returns a 2nd-order tensor. This
            %   tensor field being the result of the double contraction:
            %   lhs_{abcd}*rhs_{cd}.
            
            nElements = size(lhs, 5);

            % Memory Allocation
            newData(3,3,nElements) = rhs(1,1,1);

            ldata = lhs;
            rdata = rhs;

            for i = 1:nElements
                for a = 1:3
                    for b = 1:3
                        % Optimization needed! The product:
                        %   matrix*rdata(:,:,i) is not necessary and
                        % neither reshape.
                        matrix = reshape(ldata(a,b,:,:,i), [3, 3]);
                        newData(a,b,i) = sum(diag(matrix*rdata(:,:,i).'));
                    end
                end
            end

            newTensorField = newData;
        end
    end
end
