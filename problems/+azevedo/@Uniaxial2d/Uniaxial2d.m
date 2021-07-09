classdef Uniaxial2d < azevedo.Uniaxial & handle
    %Uniaxial2d compression/elongation
    %
    % TODO:
    %   - Elongation
    %
    % Nomenclature.: 
    % - The letters U stand for displacement vector
    % respectively.
    % - The U_o, U_oo stands for U at old, and old old time respectively.

    
    properties

    end
    
    methods
        
        function obj = Uniaxial2d(varargin)
            
            obj@azevedo.Uniaxial(varargin{:})
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
                    planeStress(true),...
                    varargin{:}...
                );
        end
        
        function obj = setBoundaryConditions(obj)
            
            for iBoundary=1:obj.mesh().numberOfBoundaries
    
                boundary = obj.mesh().boundaries(iBoundary);

                bName = boundary.userName;
                boundCondition = [];
                
                if strcmp(bName,'top') || strcmp(bName,'bottom')

                    % ion free boundary
                    boundCondition.kind = 'solidTraction';
                    %boundCondition.kind = 'solidSymmetry';

                elseif strcmp(bName,'right')

                    boundCondition.kind = 'uniformVariableDisplacement';
                    
                    boundCondition.displacementFunction = ...
                        @(time) obj.deltaLength(time);
                
                elseif strcmp(bName,'left')

                    % The base is fixed on the ground.
                    %boundCondition.kind = 'fixedDisplacement';
                    boundCondition.kind = 'solidSymmetry';
                    
                elseif strcmp(bName,'front') || strcmp(bName,'back')

                    % Set case as 2d.
                    boundCondition.kind = 'empty';
                    
                else
                    ME = MException('Uniaxial2d:boundaryCondition',...
                        'Invalid boundary name: %s', bName);
                    throw(ME)
                end
                
                boundCondition.index = iBoundary;
                obj.boundaryConditions_{iBoundary} = boundCondition;
            end
        end
        
        function setCamera(obj)
            
            xLim = 1.1;
            
            if obj.deltaLength_ > 0
                xLim = 1 + obj.deltaLength_ + 0.1;
            end
            
            set...
            (...
                gca,...
                'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],...
                'XLim',[-0.1 xLim],...
                'YLim',[-0.1 1.1],...
                'ZLim',[-0.1 1.1]...
             )
            xlabel('x') 
            ylabel('y')
            zlabel('z')
        end
    end
end

