classdef tensorField
    %TensorList It is essentially an ordered lists of tensors.
    %   Detailed explanation goes here
    
    properties
        % 3x3xN array of numbers
        data
    end
    
    methods
        
        function obj = tensorField(input)
            % tensorField An ordered lists of tensors.
            %   tensorField() Creates an uninitilised field. No memory
            %   allocation is done.
            %
            %   tensorField(N) Creates a field with N elements
            %
            %   tensorField(f) Creates a copy of f if f is a tensorField.
            %
            %   tensorField(Array) Creates a field with the content of
            %   Array ( = 3x3xN array of numbers  ). Basically it converts
            %   MATLAB's array to tensorField. 
            %   Ex.:
            %       Array = zeros(3,3,10);
            %       ... fill Array ...
            %       s = tensorField(Array);
            
            % Not an empty field
            if nargin ~= 0
              
                % Create a field of size 'input'.
                if length(input) == 1 ...
                   && (isfloat(input) || isinteger(input))
                    
                    obj.data = zeros(3, 3, input);

                % Copy constructor.
                elseif isa(input, 'tensorField')
                    
                    obj.data = input.data;

                % Copy from raw data.
                elseif size(input, 1) == 3 && size(input, 2) == 3 ...
                        && length(size(input)) == 3
                    
                    obj.data = input;
                else
                    ME = MException('tensorField:constructor',...
                        'Invalid parameters.');
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
        
        function value = plus(lhs, rhs)
            
            if isa(lhs, 'tensorField') && isa(rhs, 'tensorField')
                
                value = tensorField(lhs.data + rhs.data);
            else
                ME = MException('tensorField:plus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function value = minus(lhs, rhs)
            
            if isa(lhs, 'tensorField') && isa(rhs, 'tensorField')
                
                value = tensorField(lhs.data - rhs.data);
            else
                ME = MException('tensorField:minus',...
                        'Cannot sum elements. Not implemented yet.');
                throw(ME);
            end
        end
        
        function value = tensor_X_tensorPar(lhs,rhs)
           
            newData = ...
                reshape...
                (...
                    tensorField.mulMatxMatLinear...
                    (...
                        lhs.data(1,1,:),...
                        lhs.data(2,1,:),...
                        lhs.data(3,1,:),...
                        lhs.data(1,2,:),...
                        lhs.data(2,2,:),...
                        lhs.data(3,2,:),...
                        lhs.data(1,3,:),...
                        lhs.data(2,3,:),...
                        lhs.data(3,3,:),...
                        rhs.data(1,1,:),...
                        rhs.data(2,1,:),...
                        rhs.data(3,1,:),...
                        rhs.data(1,2,:),...
                        rhs.data(2,2,:),...
                        rhs.data(3,2,:),...
                        rhs.data(1,3,:),...
                        rhs.data(2,3,:),...
                        rhs.data(3,3,:) ...
                    ),3,3,size(lhs.data,3)...
                ); 
            
            value = tensorField(newData);
        end
        
        function value = tensor_X_tensor(lhs,rhs)
           
            newData = lhs.data;

            ldata = lhs.data;
            rdata = rhs.data;

            for i = 1:size(lhs.data,3)

                newData(:,:,i) = ldata(:,:,i) * rdata(:,:,i);
            end

            value = tensorField(newData);
        end
        
        function value = tensor_X_vectorPar(lhs,rhs)
            
            nElements = size(lhs.data,3);
            
            newData = ...
                reshape...
                (...
                    tensorField.mulMatxVecLinear...
                    (...
                        reshape(lhs.data(1,1,:),1,nElements),...
                        reshape(lhs.data(2,1,:),1,nElements),...
                        reshape(lhs.data(3,1,:),1,nElements),...
                        reshape(lhs.data(1,2,:),1,nElements),...
                        reshape(lhs.data(2,2,:),1,nElements),...
                        reshape(lhs.data(3,2,:),1,nElements),...
                        reshape(lhs.data(1,3,:),1,nElements),...
                        reshape(lhs.data(2,3,:),1,nElements),...
                        reshape(lhs.data(3,3,:),1,nElements),...
                        rhs.data(1,:),...
                        rhs.data(2,:),...
                        rhs.data(3,:) ...
                    ),3,nElements...
                );  
            
            value = vectorField(newData);
        end
        
        function value = tensor_X_vector(lhs,rhs)
             
            newData = rhs.data;

            ldata = lhs.data;
            rdata = rhs.data;

            for i = 1:size(lhs.data,3)

                newData(:,i) = ldata(:,:,i) * rdata(:,i);
            end

            value = vectorField(newData);
        end
        
        function value = tensor_X_scalarPar(lhs,rhs)
            
            nElements = size(lhs.data,3);
            
            newData = ...
                reshape...
                (...
                    tensorField.mulMatxScaLinear...
                    (...
                        lhs.data(1,1,:),...
                        lhs.data(2,1,:),...
                        lhs.data(3,1,:),...
                        lhs.data(1,2,:),...
                        lhs.data(2,2,:),...
                        lhs.data(3,2,:),...
                        lhs.data(1,3,:),...
                        lhs.data(2,3,:),...
                        lhs.data(3,3,:),...
                        reshape(rhs.data(:),1,1,nElements)...
                    ),3,3,nElements...
                );  
            
            value = tensorField(newData);
        end
        
        function value = tensor_Div_scalarPar(lhs,rhs)
            
            nElements = size(lhs.data,3);
            
            newData = ...
                reshape...
                (...
                    tensorField.divMatxScaLinear...
                    (...
                        lhs.data(1,1,:),...
                        lhs.data(2,1,:),...
                        lhs.data(3,1,:),...
                        lhs.data(1,2,:),...
                        lhs.data(2,2,:),...
                        lhs.data(3,2,:),...
                        lhs.data(1,3,:),...
                        lhs.data(2,3,:),...
                        lhs.data(3,3,:),...
                        reshape(rhs.data(:),1,1,nElements)...
                    ),3,3,nElements...
                );  
            
            value = tensorField(newData);
        end
        
        function value = mtimes(lhs, rhs)
            
            global VECTORIZE;
            
            if length(rhs) == 1 && isfloat(rhs)
                
                value = tensorField(lhs.data * rhs);
                
            elseif length(lhs) == 1 && isfloat(lhs)
                
                value = tensorField(lhs * rhs.data);
                
            elseif isa(lhs, 'tensorField') && isa(rhs, 'tensorField')
                
                if VECTORIZE
                    value = tensor_X_tensorPar(lhs,rhs);
                else
                    value = tensor_X_tensor(lhs,rhs);
                end
                
                
            elseif isa(lhs, 'tensorField') && isa(rhs, 'vectorField')
                
                if VECTORIZE
                    value = tensor_X_vectorPar(lhs,rhs);
                else
                    value = tensor_X_vector(lhs,rhs);
                end
                
            elseif isa(lhs, 'tensorField') && isa(rhs, 'scalarField')
                
                if VECTORIZE
                    value = tensor_X_scalarPar(lhs,rhs);
                else
                    newData = lhs.data;

                    ldata = lhs.data;
                    rdata = rhs.data;

                    for i = 1:size(lhs.data,3)

                        newData(:,:,i) = ldata(:,:,i) * rdata(i);
                    end

                    value = tensorField(newData);
                end
            
            elseif isa(lhs, 'scalarField') && isa(rhs, 'tensorField')
                
                newData = rhs.data;
                
                ldata = lhs.data;
                rdata = rhs.data;
                
                for i = 1:size(rdata,3)
                    
                    newData(:,:,i) = ldata(i)*rdata(:,:,i);
                end
                
                value = tensorField(newData);
                
            elseif isa(rhs, 'tensor4thOrderField')
                
                nElements = rhs.numberOfElements();
                
                % Memory allocation.
                newData = rhs.data;
                
                % For performance.
                ldata = lhs.data;
                rdata = rhs.data;
                
                for i=1:nElements
                    for j=1:3
                        newData(j,:,:,:,i) = ...
                              ldata(j,1,i)*rdata(1,:,:,:,i)...
                            + ldata(j,2,i)*rdata(2,:,:,:,i)...
                            + ldata(j,3,i)*rdata(3,:,:,:,i);
                    end
                end
                
                value = tensor4thOrderField(newData);
            else
                ME = MException('tensorField:mtimes',...
                        'Cannot mutiply elements');
                throw(ME);
            end
        end
        
        function value = mrdivide(lhs, rhs)
            
            global VECTORIZE;
            
            if length(rhs) == 1 && isfloat(rhs)
                
                value = tensorField(lhs.data / rhs);
                
            elseif isa(lhs, 'tensorField') && isa(rhs, 'scalarField')
                
                if VECTORIZE
                    value = tensor_Div_scalarPar(lhs,rhs);
                else
                    newData = lhs.data;
                    rdata = rhs.data;

                    for i = 1:size(newData,3)

                        newData(:,:,i) = newData(:,:,i) / rdata(i);
                    end

                    value = tensorField(newData);
                end
            else
                ME = MException('tensorField:mrdivide',...
                        'Cannot divide elements');
                throw(ME);
            end
        end
        
        function value = transpose(obj)
            
            value = tensorField(permute(obj.data, [2 1 3]));
        end
        
        function value = inv(obj)
            
            global VECTORIZE;
            
            if VECTORIZE
                % Calculate and convert back to matrix
                newData = ...
                    reshape...
                    (...
                        tensorField.invLinear...
                        (...
                            obj.data(1,1,:),...
                            obj.data(2,1,:),...
                            obj.data(3,1,:),...
                            obj.data(1,2,:),...
                            obj.data(2,2,:),...
                            obj.data(3,2,:),...
                            obj.data(1,3,:),...
                            obj.data(2,3,:),...
                            obj.data(3,3,:) ...
                        ),3,3,size(obj.data,3)...
                    );
            else
                newData = obj.data;
                
                for i = 1:size(newData,3)

                    newData(:,:,i) = inv(newData(:,:,i));
                end
            end
            value = tensorField(newData);
        end
        
        function value = det(obj)
            
            global VECTORIZE
            
            if VECTORIZE
                newData = ...
                    reshape(...
                        tensorField.detLinear...
                        (...
                            obj.data(1,1,:),...
                            obj.data(2,1,:),...
                            obj.data(3,1,:),...
                            obj.data(1,2,:),...
                            obj.data(2,2,:),...
                            obj.data(3,2,:),...
                            obj.data(1,3,:),...
                            obj.data(2,3,:),...
                            obj.data(3,3,:) ...
                        ),1,size(obj.data,3)...
                    );
                value = scalarField(newData);
            else
                rawData = obj.data;
                N = size(obj.data,3);
                newData = zeros(1,N);

                for i = 1:N

                    newData(i) = det(rawData(:,:,i));
                end
            
                value = scalarField(newData);
            end
        end
        
        function value = trace(obj)
            
            rawData = obj.data;
            N = size(obj.data,3);
            newData = zeros(1,N);
                
            for i = 1:N

                newData(i) = trace(rawData(:,:,i));
            end
            
            value = scalarField(newData);
        end
        
        function endVal = end(obj, k, n)
            
            if n == 1
                endVal = size(obj.data, 3);
            else
                ME = MException('tensorField:end',...
                        'Not implemented.');
                throw(ME);
                %endVal = builtin('end', obj.data, k, n);
            end
        end
        
        function nElements = numberOfElements(obj)
            
            nElements = size(obj.data, 3);
        end
        
        function output = internalField(obj, i, j, k)

            if nargin == 2

                output = obj.data(:,:,i);

            elseif nargin == 3

                output = obj.data(:,i,j);
                
            elseif nargin == 4

                output = obj.data(i,j,k);

            else
                error('Use one to three parameters.');
            end
        end
        
        function output = range(obj, range)

            if nargin == 2

                output = obj.data(:,:,range);
            else
                error('Use one parameter.');
            end
        end
        
        function output = subField(obj, range)

            output = tensorField(obj.range(range));
        end
        
        function output = field(obj, i, j, k)

            if nargin == 2

                output = obj.data(:,:,i);

            elseif nargin == 3

                output = obj.data(:,i,j);
                
            elseif nargin == 4

                output = obj.data(i,i,k);

            else
                error('Use one to three parameters.');
            end
        end
        
        function [eigVectors, eigValues] = eigen(obj)
            % EIGEN Eigenvalues and eigenvectors.
            %   [V D] = aTensorField.eigen() produces 2 tensorFields. 
            %   Each D(i) is diagonal matrix of eigenvalues and each V(i) a 
            %   full matrix whose columns are the corresponding 
            %   eigenvectors of aTensorField(i), so that 
            %   aTensorField(i)*V(i) = V(i)*D or aTensorField*V = V*D.
            
            rawVectors = obj.data;
            rawValues = obj.data;
            rawData = obj.data;
            
            for i = 1:size(rawData,3)
                
                [ rawVectors(:,:,i), rawValues(:,:,i)] = ...
                    eig(rawData(:,:,i));
            end
            
            eigValues = tensorField(rawValues);
            eigVectors = tensorField(rawVectors);
        end
        
        function row = rowAsVector(obj, idx)
            %rowAsVector Return a row as vectorField.
            %   rowAsVector(idx) returns the idx-th row as a vectorField.
            
            format = [3 obj.numberOfElements()];
            array = obj.data(idx,:,:);
            row = vectorField(reshape(array, format));
        end
        
        function col = colAsVector(obj, idx)
            %colAsVector Return a column as a vectorField.
            %   colAsVector(idx) returns the idx-th column as a
            %   vectorField.
            
            format = [3 obj.numberOfElements()];
            array = obj.data(:,idx,:);
            col = vectorField(reshape(array, format));
        end
        
        function identity = identity(obj)
            
            identity = tensorField.identity_(obj.numberOfElements());
        end
        
        function obj = zero(obj)
            
           obj.data(:,:,:) = 0;
        end
    end
    
    methods (Static)

        function identity = identity_(nElements)
            
            persistent hashMap;
            
            if isempty(hashMap)
               hashMap = ...
                   containers.Map('KeyType', 'single', 'ValueType', 'any');
            end
            
            if isKey(hashMap, nElements)
                identity = hashMap(nElements);
            else
                Idata = repmat(eye(3,3),1,1,nElements);
                identity = tensorField(Idata);
                hashMap(nElements) = identity;
            end
        end
        
        function output = ...
                invLinear(a1_1,a2_1,a3_1,a1_2,a2_2,a3_2,a1_3,a2_3,a3_3)
            % invLinear Function for vectorized inversion
            
            d = ...
                (...
                    a1_1.*a2_2.*a3_3 ...
                   -a1_1.*a2_3.*a3_2 ...
                   -a1_2.*a2_1.*a3_3 ...
                   +a1_2.*a2_3.*a3_1 ...
                   +a1_3.*a2_1.*a3_2 ...
                   -a1_3.*a2_2.*a3_1 ...
                );
            
            output = ...
                [...
                    (a2_2.*a3_3-a2_3.*a3_2)./d;...
                   -(a2_1.*a3_3-a2_3.*a3_1)./d;...
                    (a2_1.*a3_2-a2_2.*a3_1)./d;...
                   -(a1_2.*a3_3-a1_3.*a3_2)./d;...
                    (a1_1.*a3_3-a1_3.*a3_1)./d;...
                   -(a1_1.*a3_2-a1_2.*a3_1)./d;...
                    (a1_2.*a2_3-a1_3.*a2_2)./d;...
                   -(a1_1.*a2_3-a1_3.*a2_1)./d;...
                    (a1_1.*a2_2-a1_2.*a2_1)./d...
                ];
        end
        
        function output = ...
                detLinear(A1_1,A2_1,A3_1,A1_2,A2_2,A3_2,A1_3,A2_3,A3_3)
            
            output = ...
                 A1_1.*A2_2.*A3_3 ...
                -A1_1.*A2_3.*A3_2 ...
                -A1_2.*A2_1.*A3_3 ...
                +A1_2.*A2_3.*A3_1 ...
                +A1_3.*A2_1.*A3_2 ...
                -A1_3.*A2_2.*A3_1;
        end
        
        function output = ...
                mulMatxMatLinear...
                (...
                    a1_1,a2_1,a3_1,a1_2,a2_2,a3_2,a1_3,a2_3,a3_3,...
                    b1_1,b2_1,b3_1,b1_2,b2_2,b3_2,b1_3,b2_3,b3_3 ...
                )
            
            output = ...
                [...
                    a1_1.*b1_1+a1_2.*b2_1+a1_3.*b3_1;...
                    a2_1.*b1_1+a2_2.*b2_1+a2_3.*b3_1;...
                    a3_1.*b1_1+a3_2.*b2_1+a3_3.*b3_1;...
                    a1_1.*b1_2+a1_2.*b2_2+a1_3.*b3_2;...
                    a2_1.*b1_2+a2_2.*b2_2+a2_3.*b3_2;...
                    a3_1.*b1_2+a3_2.*b2_2+a3_3.*b3_2;...
                    a1_1.*b1_3+a1_2.*b2_3+a1_3.*b3_3;...
                    a2_1.*b1_3+a2_2.*b2_3+a2_3.*b3_3;...
                    a3_1.*b1_3+a3_2.*b2_3+a3_3.*b3_3 ...
                ];
        end
        
        function output = ...
                mulMatxVecLinear...
                (...
                    A1_1,A2_1,A3_1,A1_2,A2_2,A3_2,A1_3,A2_3,A3_3,...
                    b1,b2,b3...
                )
            
            output = ...
                [...
                    A1_1.*b1+A1_2.*b2+A1_3.*b3;...
                    A2_1.*b1+A2_2.*b2+A2_3.*b3;...
                    A3_1.*b1+A3_2.*b2+A3_3.*b3 ...
                ];
        end
        
        function output = ...
                mulMatxScaLinear...
                (...
                    A1_1,A2_1,A3_1,A1_2,A2_2,A3_2,A1_3,A2_3,A3_3,...
                    a...
                )
            
            output = ...
                [...
                    A1_1.*a; A2_1.*a; A3_1.*a; ...
                    A1_2.*a; A2_2.*a; A3_2.*a; ...
                    A1_3.*a; A2_3.*a; A3_3.*a; ...
                ];
        end
        
        function output = ...
                divMatxScaLinear...
                (...
                    A1_1,A2_1,A3_1,A1_2,A2_2,A3_2,A1_3,A2_3,A3_3,...
                    a...
                )
            
            output = ...
                [...
                    A1_1./a; A2_1./a; A3_1./a; ...
                    A1_2./a; A2_2./a; A3_2./a; ...
                    A1_3./a; A2_3./a; A3_3./a; ...
                ];
        end
    end
end
