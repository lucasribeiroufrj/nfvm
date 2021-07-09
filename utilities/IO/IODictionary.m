classdef IODictionary
    %IODictionary Summary of this class goes here
    %   Detailed explanation goes here
    %   ATTENTION: THIS CLASS IS A DRAFT
    %
    % Dependency:
    % - TokenStream
    % - containers.Map (MATLAB)
    %
    % /-----------------------------\
    % | Grammar for this dictionary |
    % \-----------------------------/
    %
    % item:
    % 
    %     key { item } item *
    %     key value ; item*
    % 
    % key:
    %     identifier
    % 
    % value:
    %     vector
    %     identifier
    %     nonuniform
    %     uniform
    % 
    % vector:
    %     [ number+ ]
    % 
    % nonuniform:
    %     List < vectorType > size ( vectorTuple+ )
    % 
    % uniform:
    %     ( number+ ) | number
    % 
    % vectorType:
    %     'vector'
    % 
    % vectorTuple:
    %     ( number number number )
    % 
    % size:
    %     integer
    
    properties
        
        dictionary_
    end
    
    methods
        function obj = IODictionary(input)
            %IODictionary Construct an instance of this class
            %   Detailed explanation goes here
            %
            %   IODictionary(tokenStream) is a IODictionary created using
            %   an instance of the class TokenStream.
            %
            %   IODictionary(fileDirectory) is a IODictionary created using
            %   using a string with the name of a file.
            
            if isa(input, 'TokenStream')
                
                tokenStream = input;
            else
                fileDirectory = input;
                tokenStream = TokenStream(fileDirectory, true);
            end
            
            obj.dictionary_ = IODictionary.createDictionary(tokenStream);
        end
        
        function disp(obj)
            
            IODictionary.print(obj.dictionary_);
        end
        
        function map = getMap(obj)
            
            map = obj.dictionary_;
        end
        
        function varargout = subsref(obj, S)
            
            switch S(1).type
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                    
                case '()'
                    if length(S) == 1
                        
                        if length(S.subs) == 1
                            
                            [varargout{1:nargout}] = ...
                                obj.dictionary_(S.subs{1});
                        else
                            error('Not a valid indexing expression')
                        end
                    end
                otherwise
                    error('Not a valid indexing expression')
            end
        end
    end
    
    methods (Static)
       
        function dictionary = createDictionary(tokenStream)
            %createDictionary Summary of this method goes here
            %   Detailed explanation goes here
            
            dictionary = [];

            while true

                token = tokenStream.get();

                if isempty(token)
                    return;
                end

                if token == '}'
                    tokenStream.putTokenBack(token);
                    return
                end

                dictionary = IODictionary.createMap();
                nextToken = tokenStream.get();

                if nextToken == '{'
                    
                    tokenStream.putTokenBack(token);
                    key = IODictionary.key(tokenStream);
                    
                    subDictionary = ...
                        IODictionary.createDictionary(tokenStream);
                    dictionary(key) = subDictionary;

                    IODictionary.consumeRightCurlyBracket(tokenStream);
                else

                    tokenStream.putTokenBack(nextToken);
                    value = IODictionary.value(tokenStream);
                    
                    tokenStream.putTokenBack(token);
                    key = IODictionary.key(tokenStream);
                    
                    dictionary(key) = value;

                    IODictionary.consumeSemicolon(tokenStream);
                end

                item = IODictionary.createDictionary(tokenStream);

                if isempty(item)
                    return
                end

                dictionary = [ dictionary; item ]; %#ok<AGROW>
            end
        end
        
        function map = createMap()
            %createMap Summary of this method goes here
            %   Detailed explanation goes here
            
            map = containers.Map('KeyType','char','ValueType','any');
        end

        function key = key(tokenStream)
            %key Summary of this method goes here
            %   Detailed explanation goes here

            identifier = IODictionary.identifier(tokenStream);
            key = identifier;
        end

        function identifier = identifier(tokenStream)
            %identifier Summary of this method goes here
            %   Detailed explanation goes here

            identifier = tokenStream.get();
        end

        function value = value(tokenStream)
            %value Summary of this method goes here
            %   Detailed explanation goes here

            token = tokenStream.get();

            if token == '['

                value = IODictionary.vector(tokenStream);

            elseif strcmp(token, 'nonuniform')

                value = IODictionary.nonuniform(tokenStream);

            elseif strcmp(token, 'uniform')

                value = IODictionary.uniform(tokenStream);
            else
                tokenStream.putTokenBack(token);
                value = IODictionary.identifier(tokenStream);
            end
        end

        function vector = vector(tokenStream)
            %vector Summary of this method goes here
            %   Detailed explanation goes here

            number = tokenStream.get();
            vector = number;

            while true
                token = tokenStream.get();

                if token == ']'
                    break;
                end

                number = token;
                vector = [ vector number ]; %#ok<AGROW>
            end
        end
        
        function consumeChar(tokenStream, anyChar)
            
            token = tokenStream.get();

            if token ~= anyChar

                error('Missing %s ', anyChar);
            end
        end

        function consumeRightCurlyBracket(tokenStream)
            %consumeRightCurlyBracket Summary of this method goes here
            %   Detailed explanation goes here

            IODictionary.consumeChar(tokenStream, '}');
        end

        function consumeSemicolon(tokenStream)
            %consumeSemicolon Summary of this method goes here
            %   Detailed explanation goes here

            IODictionary.consumeChar(tokenStream, ';');
        end

        function consumeLeftParenthesis(tokenStream)
            %consumeLeftParenthesis Summary of this method goes here
            %   Detailed explanation goes here

            IODictionary.consumeChar(tokenStream, '(');
        end

        function consumeRightParenthesis(tokenStream)
            %consumeRightParenthesis Summary of this method goes here
            %   Detailed explanation goes here

            IODictionary.consumeChar(tokenStream, ')');
        end

        function nonuniform = nonuniform(tokenStream)
            %nonuniform Summary of this method goes here
            %   Detailed explanation goes here

            token = tokenStream.get();

            if strcmp(token, 'List')

                token = tokenStream.get();

                if token ~= '<'
                    error('Missing token <')
                end

                token = tokenStream.get();

                if ~strcmp(token, 'vector')

                    error('vector is the only token supported.');
                end

                token = tokenStream.get();

                if token ~= '>'
                    error('Missing token <')
                end

                token = tokenStream.get();
                listSize = token;
                tuples = zeros(3,listSize);

                IODictionary.consumeLeftParenthesis(tokenStream);

                for i = 1:listSize

                    tuple = IODictionary.vectorTuple(tokenStream);
                    tuples(:,i) = tuple;
                end

                IODictionary.consumeRightParenthesis(tokenStream);
            else
                error('List is the only token supported.');
            end

            nonuniform = tuples;
        end

        function uniform = uniform(tokenStream)
            %uniform Summary of this method goes here
            %   Detailed explanation goes here
            
            token = tokenStream.get();

            if token == '('

                %IODictionary.consumeLeftParenthesis(tokenStream);
                uniform = [];

                while true
                    token = tokenStream.get();

                    if token == ')'
                        break;
                    end

                    tokenStream.putTokenBack(token);
                    number = IODictionary.number(tokenStream);

                    uniform = [ uniform number ]; %#ok<AGROW>
                end
            else
                tokenStream.putTokenBack(token);
                uniform = IODictionary.number(tokenStream);
            end
        end

        function vectorTuple = vectorTuple(tokenStream)
            %vectorTuple Summary of this method goes here
            %   Detailed explanation goes here

            IODictionary.consumeLeftParenthesis(tokenStream);

            vectorTuple = [ ...
                    IODictionary.number(tokenStream), ...
                    IODictionary.number(tokenStream), ...
                    IODictionary.number(tokenStream)  ...
                ];

            IODictionary.consumeRightParenthesis(tokenStream);
        end

        function number = number(tokenStream)
            %number Summary of this method goes here
            %   Detailed explanation goes here

            token = tokenStream.get();

            if token == '-'
                token = -tokenStream.get();
            end

            number = token;
        end
        
        function print(dictionary, tabs)
            %print Summary of this method goes here
            %   Detailed explanation goes here

            if nargin == 1
                tabs = '';
            end
            
            keys = dictionary.keys;
            values = dictionary.values;

            for i = 1:length(dictionary)

                IODictionary.printItem(keys(i), values(i), tabs);
            end
        end

        function printItem(key, value, tabs)
            %printItem Summary of this method goes here
            %   Detailed explanation goes here

            fprintf('%s%s', tabs, key{1})

            if isa(value{1},'containers.Map')

                fprintf('\n%s{\n', tabs)

                tabs = [ tabs, '    '];
                IODictionary.print(value{1}, tabs);

                fprintf('%s}\n', tabs(1:length(tabs)-4))
            else
                fprintf('   ');
                disp(value{1});
            end
        end
    end
end

