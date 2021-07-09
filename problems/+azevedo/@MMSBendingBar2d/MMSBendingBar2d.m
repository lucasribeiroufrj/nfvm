classdef MMSBendingBar2d < azevedo.MMSBendingBar & handle
    %MMSBendingBar2d See [7, Kamojjala - 2013].
    %
    % TODO:
    %   - Traction-based boundary condition.
    %
    % Nomenclature.: 
    % - The letters U stand for displacement vector
    % respectively.
    % - The U_o, U_oo stands for U at old, and old old time respectively.
    %
    % Ref.:
    % [Kamojjala - 2013] - Verification tests in solid mechanics.
    
    properties
        
    end
    
    methods
        
        function obj = MMSBendingBar2d(varargin)
            
            obj@azevedo.MMSBendingBar(varargin{:})
        end
        
        function obj = setBoundaryConditions(obj)
            
            for iBoundary=1:obj.mesh().numberOfBoundaries
    
                boundary = obj.mesh().boundaries(iBoundary);
                bName = boundary.userName;
                
                if strcmp(bName,'front') || strcmp(bName,'back')
                    
                    % Set case as 2d.
                    kind = 'empty';
                elseif strcmp(bName,'bottom')

                    % The base is fixed on the ground.
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
        
        function material = createMaterial...
                (...
                    obj,...
                    materialEnum,...
                    varargin...
                )
            
            material = ...
                materials.createMaterial...
                (...
                    materialEnum,...
                    planeStress(false),...
                    varargin{:}...
                );
        end
        
        function setCamera(obj)
            
            margin = 0.2*obj.barHeight_;
            yLim = obj.barHeight_ + margin*3;
            xLim = 0.9*obj.barHeight_;%1.2;
            
%             if obj.barHeight_ > 0
%                 xLim = 1 + obj.deltaLength_ + 0.1;
%             end
            
%             set...
%             (...
%                 gca,...
%                 'DataAspectRatio',[1 1 1],...
%                 'PlotBoxAspectRatio',[1 1 1],...
%                 'XLim',[-obj.barHeight_*2 xLim],...
%                 'YLim',[-0.2 obj.barHeight_*2.0],...
%                 'ZLim',[-0.2 1.2]...
%              )
%             xlabel('x') 
%             ylabel('y')
%             zlabel('z')

            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],...
                'XLim',[-xLim 1 + margin],...
                'YLim',[-margin yLim],...
                'ZLim',[-margin 1 + margin]...
            )
        end
    end
end

