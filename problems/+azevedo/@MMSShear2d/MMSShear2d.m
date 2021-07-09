classdef MMSShear2d < azevedo.MMSShear & handle
    %MMSShear Shear case with analytical answer (Method of Manufactured
    %   Solutions).
    %
    % Ref.:
    % [Bonet - 2008] - Nonlinear continuum mechanics for finite element
    % analysis
    
    properties

    end
    
    methods
        
        function obj = MMSShear2d(varargin)
            
            obj@azevedo.MMSShear(varargin{:})
        end
        
        function obj = setBoundaryConditions(obj)
            
            for iBoundary=1:obj.mesh().numberOfBoundaries
    
                boundary = obj.mesh().boundaries(iBoundary);
                bName = boundary.userName;
                
                if strcmp(bName,'back') || strcmp(bName,'front')
                    
                    kind = 'empty';

                elseif strcmp(bName,'bottom')

                    %kind = 'MMSDisplacement';
                    kind = 'fixedDisplacement';
                else
                    if obj.useTraction_

                        kind = 'MMSTraction';
                    else
                        kind = 'MMSDisplacement';
                    end
                end
                
                obj.boundaryConditions_{iBoundary}.index = iBoundary;
                obj.boundaryConditions_{iBoundary}.problem = obj;
                obj.boundaryConditions_{iBoundary}.kind = kind;
            end
        end
        
        function setCamera(obj)
            
            xLim = 1.2;
            
            if obj.shearFactor_ > 0
                xLim = 1 + obj.shearFactor_ + 0.1;
            end
            
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
    
    methods(Static)
       
        function fullPaths = getMeshPaths()
            %getMeshPaths return a cell with all mesh paths.

            fullpath = mfilename('fullpath');
            className = mfilename;
            fullPaths = Problem.getMeshPaths(fullpath, className);
        end
    end
end

