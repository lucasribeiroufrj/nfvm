classdef (Abstract) AbstractScalarField
    %AbstractScalarField Interface class for a scalar field
    %   Let scalarField be any from-AbstractScalarField derived class, then
    % scalarField is ordered lists of scalars.
    %   It must have a constructor which having the following contract:
    %
    %   scalarField() Creates an uninitilised field. No memory
    %   allocation is done.
    %
    %   scalarField(N) Creates a field with N elements
    %
    %   scalarField(f) Creates a copy of f if f is a scalarField.
    %
    %   scalarField(Array) Creates a field with the content of
    %   Array ( = 1xN array of numbers  ). 
    %   Ex.:
    %       Array = zeros(1,10);
    %       ... fill Array ...
    %       s = scalarField(Array);
    
    properties (Abstract)
        %(Access = {...
        %    ?AbstractTensorField,...
        %    ?AbstractTensorField}...
        %)
        
        % An array, i.e. an ordered list of scalars
        data
    end
    
    methods (Abstract, Access=public)
        
        % Should return the length of data member
        nElements = numberOfElements(obj)
        
        % BEG - Operator overloading
        %
        % They are all element-wise operations
        value = mtimes(lhs, rhs)
        
        value = plus(lhs, rhs)
        
        value = minus(lhs, rhs)
        
        value = mrdivide(lhs, rhs)
        
        value = uminus(obj)
        
        value = mpower(obj, exponent)

        value = inv(obj)
        % END - Operator overloading
        
        
        % BEG - Accessory methods
        %
        % Should return the containing array as an array of double
        array = double(anyScalarField)
        % END - Accessory methods
    end
    
    methods(Static, Abstract)
        
        % Element-wise operation
        value = sqrt(field)
        
        % Element-wise operation
        value = log(field)
    end
end

