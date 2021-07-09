classdef PrismaticBar < Problem & handle
    %PrismaticBar Summary of this class goes here
    %   Detailed explanation goes here
    %
    % Case described in [4.2; Demirdzic - 1994].
    %
    % Ref.:
    %   [Demirdzic - 1994] - Finite volume method for stress analysis in 
    % complex domains
    
    properties
        
        % Use traction or displacamente for top boundary.
        useTraction_ = true;
        
        gravity_ = [ 0; 0; -10 ];
        
        % The bar's length.
        l_ = 2.5;
        
        % To draw to mesh
        meshDrawer_;
        
        % Boundary conditions for displacement U 
        boundaryConditions_;
    end
    
    methods
        
        function obj = PrismaticBar(casePath, params)
            %PrismaticBar Construct the case described in 
            %   [4.2; Demirdzic - 1994], i.e., prismatic bar stretched by 
            %   its own weight
            
            params.material = ...
                tryGetOrDefault(params, 'material', materials.List.Steel);
            
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                    'nTimeSteps', 0 ... % Static case
                );
            
            obj@Problem(setup, params);
            
            obj.useTraction_ = tryGetOrDefault(params,'useTraction',true);
            
            if obj.useTraction_
                
                disp('Using traction boundary for top face.');
            end
            
            obj.applyUnderRelaxation(...
                tryGetOrDefault(params, 'applyUnderRelaxation', true));
            obj.setBoundaryConditions(params.material);
            
            obj.setMeshDrawer();
            
            obj.applyUnderRelaxation(...
                tryGetOrDefault(params, 'applyUnderRelaxation', true));
            
            obj.setRelaxation(...
                tryGetOrDefault(params, 'relaxation', 0.99));
            
            obj.setAlternativeTolerance(...
                tryGetOrDefault(params, 'alternativeTolerance', 1e-6));
            
            obj.setNumberOfCorrectors(...
                tryGetOrDefault(params, 'numberOfCorrectors', 150));
            
            obj.setSolidModel(params);
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
        
        function obj = setBoundaryConditions(obj, materialEnum)
           
            material = obj.createMaterial(materialEnum);
            
            for iBoundary=1:obj.mesh().numberOfBoundaries
               
                boundary = obj.mesh().boundaries(iBoundary);

                if strcmp(boundary.type,'empty')
                    continue
                end

                bName = boundary.userName;
                boundCondition = [];
                
                if strcmp(bName,'front') 
                    
                    if obj.useTraction_
                        
                        % Reaction traction caused by the gravity force.
                        boundCondition.kind = 'solidTraction';
                        density = material.density();
                        g = obj.gravity_;
                        traction = -density*g*obj.l_;
                        data = repmat(traction, [1,boundary.numberOfBFaces]);
                        boundCondition.traction = vectorField(data);
                    else
                        % Zero displacement.
                        boundCondition.kind = 'fixedDisplacement';
                    end
                    
                elseif strcmp(bName,'bottom') || strcmp(bName,'left') ...
                    || strcmp(bName,'back') || strcmp(bName,'right') ...
                    || strcmp(bName,'top')
                
                    % Traction free boundary
                    boundCondition.kind = 'solidTraction';
                end
                
                boundCondition.index = iBoundary;
                obj.boundaryConditions_{iBoundary} = boundCondition;
            end
        end
        
        function obj = setMeshDrawer(obj)
            
            obj.meshDrawer_ = BoundaryMeshDrawer(obj.mesh());
        end
        
        function U = U(obj)
            %U Displacement field
            
            U = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function U_o = U_o(obj)
            %Uo Displacement field at time "old"
            
            U_o = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function U_oo = U_oo(obj)
            %Uoo Displacement field at time "old-old"
            
            U_oo = volVectorField(obj.mesh(), obj.boundaryConditions_);
        end
        
        function meshDrawer = meshDrawer(obj)
            
            meshDrawer = obj.meshDrawer_;
        end
        
        function setCamera(obj)
            
            set...
            (   gca,...
                'CameraPosition', [-13.5 -16.5 4.6],...
                'CameraUpVector', [0.09 0.12 0.98],...
                'CameraViewAngle', 6.608,...
                'DataAspectRatio', [1 1 1],...
                'PlotBoxAspectRatio', [1 1 1]...
            );
        end
        
        function bool = hasBodyForce(obj)
           
            bool = true;
        end
        
        function RHS = bodyForce(obj)
            
            force = obj.getBodyForce();
            
            RHS = -force.data.';
        end
        
        function force = getBodyForce(obj)
            
            mesh = obj.mesh();
            N = mesh.numberOfElements();
            
            density = obj.materialLaw().density();
            data = repmat(density*obj.gravity_, [1,N]);
            density_g = vectorField(data);
            
            force = mesh.volume() * density_g;
        end
        
        function phi = analyticalDisplacement(obj, mesh)
            
            Mu = obj.materialLaw().mu();
            Lambda = obj.materialLaw().lambda();

            [ E, Nu ] = ...
                computeLameConstants(struct('Lambda', Lambda, 'Mu', Mu));
            density = obj.materialLaw().density();
            l = obj.l_;
            g = norm(obj.gravity_);
            
            rawPhi = zeros(3,mesh.numberOfElements);
            
            for iElement = 1:mesh.numberOfElements
                
                centroid = mesh.elements(iElement).centroid;
                x = centroid(1);
                y = centroid(2);
                z = centroid(3);
                rawPhi(1,iElement) = - Nu*density*g*x*z / E;
                rawPhi(2,iElement) = - Nu*density*g*y*z / E;
                rawPhi(3,iElement) = (density*g)/(2*E) * (z*z - l*l) + ...
                    + Nu*density*g / (2*E) * ( x*x + y*y );
            end
            
            phi = vectorField(rawPhi);
        end
        
        function maxDisplacementError(obj, numericalVolU)
             
            mesh = obj.mesh();
            
            infoNumeric = [];
            infoAnalyti = [];
            
            nPhi = numericalVolU.internalField();
            aPhi = obj.analyticalDisplacement(mesh);
            
            for iElem = 1:mesh.numberOfElements
               
                numerical = nPhi.range(iElem);
                analytical = aPhi.range(iElem);
                
                element = mesh.elements(iElem);
                if norm( [0.1; -0.1] - element.centroid(1:2,1) ) < 5e-4
                    z = element.centroid(3);
                    u = numerical(1);
                    v = numerical(2);
                    w = numerical(3);
                    
                    a = analytical(1);
                    b = analytical(2);
                    c = analytical(3);
                    
                    infoNumeric = [ infoNumeric [z; u; v; w] ];
                    infoAnalyti = [ infoAnalyti [z; a; b; c] ];
                end
            end
            
            % Solution + Constant is also a solution
            infoNumeric(4,:) = infoNumeric(4,:) ...
                -max(infoNumeric(4,:)-infoAnalyti(4,end));
            
            figHandle = figure(7);
            subplot(1,2,1);
            plot(infoNumeric(1,:),infoNumeric(2,:), 's');
            
            xlabel('z (m)');
            ylabel('u,v (m)');
            
            hold on;
            plot(infoAnalyti(1,:),infoAnalyti(2,:));
            
            legend('Numerical','Analytical');
            legend('Location','southwest');
            plot(infoNumeric(1,:),infoNumeric(3,:), 's');
            plot(infoAnalyti(1,:),infoAnalyti(3,:));
            
            hold off;
            subplot(1,2,2);
            plot(infoNumeric(1,:),infoNumeric(4,:), 's');
            
            xlabel('z (m)');
            ylabel('w (m)');
            
            hold on;
            subplot(1,2,2);
            plot(infoAnalyti(1,:),infoAnalyti(4,:));
            hold off;
            legend('Numerical','Analytical');
            legend('Location','southeast');
            
            set...
            (...
                figHandle,...
                'NumberTitle', 'off', ...
                'Name', ...
                'Comparison with analytical solution [4.2; Demirdzic - 1994]'...
            );
        end
        
        function obj = statistics(obj)
            
            volU = obj.solidModel.volU();
            obj.maxDisplacementError(volU);
        end
        
        function bool = includeAcceleration(obj)
            % includeAcceleration Override this function if you need to
            %   remove the temporal term acceleration. Change the sentence
            %   below to "bool = false".
            
            bool = false;
        end
    end
end

