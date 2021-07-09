classdef TokenStream < handle
    %TOKENSTREAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        fileDirectory_
        skipHeader_
        fileID_
        tokenBuffer_
    end
    
    methods
        function obj = TokenStream(fileDirectory, skipHeader)
            
            obj.fileDirectory_ = fileDirectory;
            obj.skipHeader_ = skipHeader;
            obj.fileID_ = -1;
            obj.tokenBuffer_ = [];
        end
        
        function token = get(obj)
           
            %% File already open?
            if obj.fileID_ == -1
                
                obj.fileID_ = fopen(obj.fileDirectory_, 'r');
                
                if obj.fileID_ == -1
                    error('File not found');
                end
            end
            %%
            
            if ~isempty(obj.tokenBuffer_)
               
                token = obj.tokenBuffer_;
                obj.tokenBuffer_ = [];
                return;
            end
            
            currentPosition = ftell(obj.fileID_);
            
            %% Skip header
            if currentPosition == 0 && obj.skipHeader_
                
                for i=1:16
                    fgetl(obj.fileID_);
                end
                
                currentPosition = ftell(obj.fileID_);
            end
            %%
            
            % Get only one character to test
            [c, newPosition] = textscan(obj.fileID_, '%c', 1);
    
            if isempty(c{1})
                
                token = [];
                return
            end
            
            switch c{1}
                case {';', '(', ')', '/', '{', '}', '+', '-', '[', ']',...
                        '<','>' '*', '\'}
                    token = c{1};
                    
                    %% Is this a comment?
                    if token == '/'
                        
                        %currentPosition = ftell(obj.fileID_);
                        [cc, newPosition] = textscan(obj.fileID_, '%c', 1);
                        
                        if cc{1} == '/'
                            line = fgetl(obj.fileID_);
                            token = obj.get();
                            return;
                        end
                        
                        obj.putBack(currentPosition, newPosition);
                    end
                    
                case {'.', '0', '1', '2', '3', '4', '5', '6', '7', '8',...
                        '9'}
                    % Put back the char
                    obj.putBack(currentPosition, newPosition);
                    number = textscan(obj.fileID_, '%f',1);
                    token = number{1};
                    
                case '"'
                    str = [];
                    c = textscan(obj.fileID_, '%c', 1);
                    
                    while c{1} ~= '"'
                        
                        str = [ str c{1} ]; %#ok<AGROW>
                        
                        c = textscan(obj.fileID_, '%c', 1);
                    end
                    
                    token = str;
                    
                otherwise
                    % Put back the char
                    obj.move(-1);
                    currentPosition = ftell(obj.fileID_);
                    [str, newPosition] = textscan(obj.fileID_, '%s', 1);
                    obj.putBack(currentPosition, newPosition);
                    identifier = regexp(str{1}, '[a-z_A-Z]+', ...
                        'match', 'once');
                    offset = length(identifier{1});
                    obj.move(offset);
                    token = identifier{1};
            end
        end
        
        function obj = putBack(obj, currentPosition, newPosition)
            
            offset = newPosition - currentPosition;
            fseek(obj.fileID_, -offset, 0); 
        end
        
        function obj = move(obj, offset)
            
            fseek(obj.fileID_, offset, 0); 
        end
        
        function obj = putTokenBack(obj, token)
           
            obj.tokenBuffer_ = token;
        end
        
        function delete(obj)
            
            if obj.fileID_ ~= -1
                
                fclose(obj.fileID_);
            end
        end
    end
    
end

