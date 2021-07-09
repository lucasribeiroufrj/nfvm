classdef Uniaxial2d < kamojjala.Uniaxial & handle
    %Uniaxial2d Uniaxial compression/elongation in 2d
    %   Case described in [5, Kamojjala - 2013].
    %
    % TODO:
    %   -/-
    %
    % Nomenclature.: 
    %   -/-
    %
    % Ref.:
    %   [Kamojjala - 2013] - Verification tests in solid mechanics
    
    properties
        
        emptyBoundariesId_;
        
        noSymmetryBoundaries_;
    end
    
    methods
        
        function obj = Uniaxial2d(varargin)
            
            obj@kamojjala.Uniaxial(varargin{:});
        end
        
        function obj = setBoundaryConditions(obj)
            
            mesh = obj.mesh();
            
            for iBoundary = 1:obj.mesh().numberOfBoundaries()
                
                boundary = mesh.boundaries(iBoundary);
                bName = boundary.userName;
                
                if strcmp(bName,'back') || strcmp(bName,'front')
                    
                    % Set case as 2d.
                    kind = 'empty';

                elseif strcmp(bName,'left')

                    kind = 'MMSDisplacement';
                    %kind = 'solidSymmetry';
                    
                elseif (strcmp(bName,'top') || strcmp(bName,'bottom')) ...
                        && obj.useSymmetryBoundaries_
                    
                    kind = 'solidSymmetry';
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
        
%         function gradU = analyticalVolGradU(obj, time)
%             %analyticalVolGradU analytical gradient gradU.
%             %   analyticalVolGradU(time) returns gradU at a specific
%             %   time as a volTensorField.
%             
%             gradU = analyticalVolGradU@kamojjala.Uniaxial(obj,time);
%             
%             for iBoundary = obj.emptyBoundariesId_
%                 
%                 gradU.boundaryField(iBoundary) = ...
%                     gradU.boundaryField(iBoundary).zero;
%             end
%         end
        
        function setCamera(obj)
            
            dx = 0;
            if obj.Lambda_> 1
                dx = obj.Lambda_ - 1;
            end
            
            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],...
                'XLim',[-0.1 1.1 + dx],...
                'YLim',[-0.1 1.1],...
                'ZLim',[-0.1 1.1]...
             )
            xlabel('x') 
            ylabel('y')
            zlabel('z')
        end
    end
    
    methods (Static)
        
    end
end

