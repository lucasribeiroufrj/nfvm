classdef MMSSpring2d < azevedo.MMSSpring & handle
    %MMSSpring Case similar to [imechanica] with analytical answer 
    %   (Method of Manufactured Solutions). 
    %
    % TODO:
    %   - 
    %
    % Nomenclature.: 
    % - The letters U stand for displacement vector
    % respectively.
    % - The U_o, U_oo stands for U at old, and old old time respectively.
    %
    % Ref.:
    % [imechanica] - https://imechanica.org/node/1357
    
    properties

    end
    
    methods
        
        function obj = MMSSpring2d(varargin)
            
            obj@azevedo.MMSSpring(varargin{:})
        end
        
        function obj = setBoundaryConditions(obj)
            
            for iBoundary=1:obj.mesh().numberOfBoundaries
                
                boundary = obj.mesh().boundaries(iBoundary);
                bName = boundary.userName;
                boundCondition = [];
                
                if strcmp(bName,'front') || strcmp(bName,'back')
                    
                    % Set case as 2d.
                    boundCondition.kind = 'empty';
                    
                elseif strcmp(bName,'left') || strcmp(bName,'right') ...
                        || strcmp(bName,'top') || strcmp(bName,'bottom')

                    boundCondition.kind = 'MMSDisplacement';
                    boundCondition.problem = obj;
                else
                    ME = MException('MMSSpring:boundaryCondition',...
                        'Invalid boundary name: %s', bName);
                    throw(ME)
                end
                
                boundCondition.index = iBoundary;
                obj.boundaryConditions_{iBoundary} = boundCondition;
            end
        end
        
        function setCamera(obj)
            
            xLim = 2.2;
            
            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],...
                'XLim',[-0.2 xLim],...
                'YLim',[-0.2 1.2],...
                'ZLim',[-0.2 1.2]...
             )
            xlabel('x') 
            ylabel('y')
            zlabel('z')
        end
    end
end

