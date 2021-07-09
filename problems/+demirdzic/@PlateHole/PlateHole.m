classdef PlateHole < Problem & handle
    %PlateHole Summary of this class goes here
    %   Detailed explanation goes here
    %
    % Case described in [4.3; Demirdzic - 1994].
    %
    % Ref.:
    %   [Demirdzic - 1994] - Finite volume method for stress analysis in 
    % complex domains
    
    properties
        
        % To draw to mesh.
        meshDrawer_;
        
        % Boundary conditions for displacement U.
        boundaryConditions_; 
        
        % Analytical stress as volume stress field.
        analyticalVolStressField_,
        
        % Number of times statistic function has been called.
        statisticsCounter_ = 0;
        
        % Controls when stress contour should be ploted.
        statisticFrequency_ = 5;
    end
    
    methods
        
        function obj = PlateHole(casePath, params)
            
            setup = ...
                struct...
                (...
                    'casePath', casePath,...
                    'nTimeSteps', 0 ... % Static case
                );
            
            obj@Problem(setup, params);
            
            obj.setBoundaryConditions();
            
            obj.setMeshDrawer();
            
            obj.applyUnderRelaxation(...
                tryGetOrDefault(params, 'applyUnderRelaxation', false));
            
            obj.setRelaxation(...
                tryGetOrDefault(params, 'relaxation', 1.0));
            
            obj.setNumberOfCorrectors(...
                tryGetOrDefault(params, 'numberOfCorrectors', 100000));
            
            obj.setSolutionTolerance(...
                tryGetOrDefault(params, 'solutionTolerance', 1e-6));
            
            obj.setAlternativeTolerance(...
                tryGetOrDefault(params, 'alternativeTolerance', 5e-7));
            
            obj.setSolidModel(params);
        end
        
        %% BEG - Data member access        
        function obj = setBoundaryConditions(obj)
            
            mesh = obj.mesh();
            
            obj.analyticalVolStressField_ = ...
                demirdzic.PlateHole.volStressField(obj.mesh());
            
            for iBoundary=1:mesh.numberOfBoundaries
               
                boundary = mesh.boundaries(iBoundary);
                
                bName = boundary.userName;
                boundCondition = [];
                
                if strcmp(bName,'bottom') || strcmp(bName,'left')
                    
                    boundCondition.kind = 'solidSymmetry';
                    
                elseif strcmp(bName,'top') || strcmp(bName,'right')
                    
                    boundCondition.kind = 'solidTraction';
                    boundCondition.traction = ...
                        obj.tractionAtBoundary(iBoundary);
                    
                elseif strcmp(bName,'hole')
                    
                    % Traction free boundary (traction = 0 is the default).
                    boundCondition.kind = 'solidTraction';
                    %boundCondition.traction = 0;
                    
                elseif strcmp(bName,'frontAndBack')
                    
                    boundCondition.kind = 'empty';
                else
                    ME = MException('PlateHole:setBoundaryCondition',...
                        'Invalid boundary name: %s', bName);
                    throw(ME)
                end
                
                boundCondition.index = iBoundary;
                obj.boundaryConditions_{iBoundary} = boundCondition;
            end
        end
        
        function meshDrawer = meshDrawer(obj)
            
            meshDrawer = obj.meshDrawer_;
        end
        
        function obj = setMeshDrawer(obj)
            
            obj.meshDrawer_ = BoundaryMeshDrawer(obj.mesh());
        end
        %% END
        
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
        
        function maxStressError(obj, numericalVolStressField)
            
            mesh = obj.mesh();
            
            nVolStress = numericalVolStressField.internalField();
            aVolStress = obj.analyticalVolStressField_.internalField();
            
            maxStresses = zeros(1,3);
            maxAnalytical = zeros(1,3);
            
            for iElem = 1:mesh.numberOfElements
                
                numerical = nVolStress.range(iElem);
                analytical = aVolStress.range(iElem);
                diffStress = abs(numerical - analytical);
                
                vec = [ diffStress(1,1) diffStress(1,2) diffStress(2,2) ];
                sta = [ analytical(1,1) analytical(1,2) analytical(2,2) ];
                
                for i = 1:3
                    
                    if vec(i) > maxStresses(i)
                        
                        maxStresses(i) = vec(i);
                        maxAnalytical(i) = abs(sta(i));
                    end
                end
            end
            
            sxx = sprintf('xx = %0.2e (%0.2g%%)', vec(1), ...
                maxStresses(1)/maxAnalytical(1)*100);
            sxy = sprintf('xy = %0.2e (%0.2g%%)', vec(2), ...
                maxStresses(2)/maxAnalytical(2)*100);
            syy = sprintf('yy = %0.2e (%0.2g%%)', vec(3), ...
                maxStresses(3)/maxAnalytical(3)*100);
            
            fprintf('Max stress error (Pa):\n%s\n%s\n%s\n', sxx, sxy, syy);
        end
        
        function tractionField = tractionAtBoundary(obj, iBoundary)
            %tractionAtBoundary Analytical traction. Primarely used to set 
            %   boundary condition.
            %   tractionAtBoundary(mesh, iBoundary) returns the analytical
            %   traction as vectorField over the mesh.boundaries(iBoundary)
            %   boundary.
            
            mesh = obj.mesh();
            vStress = obj.analyticalVolStressField_;
            stress = vStress.boundaryField(iBoundary);
            nf = mesh.nf().boundaryField(iBoundary);
            tractionField = vectorField(stress * nf);
        end
        
        function obj = statistics(obj)
            
            gradU = obj.solidModel.volGradU();
            mechanicalModel = obj.solidModel.mechanicalModel();
            % get the second-piola stress
            
            Piola = mechanicalModel.computePiolaField(gradU);
            
            obj.maxStressError(Piola);
            
            nDivisions = 30;
            
            if mod(obj.statisticsCounter_,obj.statisticFrequency_) == 0
                
                demirdzic.PlateHole.plotContours(Piola, nDivisions);
                obj.statisticsCounter_ = 1;
                
                U = obj.solidModel.volU();
                figure(4);
                plotDisplacement...
                    (...
                        obj.mesh(),...
                        U.internalField().data(1:2,1:end),...
                        18 ...
                    );
            else
               obj.statisticsCounter_ = obj.statisticsCounter_ + 1; 
            end
        end
        
        function bool = includeAcceleration(obj)
            % includeAcceleration Override this function if you need to
            %   remove the temporal term acceleration. Change the sentence
            %   below to "bool = false".
            
            bool = false;
        end
    end
    
    methods (Static)
        
        function faceStress = faceStressField(mesh, iBoundary)
            %faceStressField Stress at the faces as a tensorField.
            %   faceStressField(mesh, iBoundary) returns a tensorField over
            %   the mesh.boundaries(iBoundary) boundary.
        
            boundary = mesh.boundaries(iBoundary);
            faceStress = tensorField(boundary.numberOfBFaces);
            
            startFace = boundary.startFace;
            endFace = startFace + boundary.numberOfBFaces - 1;

            for iFace=startFace:endFace

                face = mesh.faces(iFace);
                idx = iFace - startFace + 1;
                faceStress(idx) = demirdzic.PlateHole.stressElem(face);
            end
        end
        
        function value = stressElem(elem)
            
            x = elem.centroid(1);
            y = elem.centroid(2);
            stress = demirdzic.PlateHole.stressHandle(x,y);
            value = [ stress [0; 0]; [0 0 0]];
        end
        
        function stress = stressHandle(x,y)

            fx = 1e+4;
            r = (x.^2 + y.^2).^(1/2);
            theta = atan(y./x);
            Q = (0.5/r).^2;
            R = 3/2*(Q.^2) .* cos(4*theta);
            S = 3/2*(Q.^2) .* sin(4*theta);
            term1 = Q .* (3/2*cos(2*theta) + cos(4*theta));
            term2 = Q .* (1/2*cos(2*theta) - cos(4*theta));
            term3 = Q .* (1/2*sin(2*theta) + sin(4*theta));
            stressXX = fx * (1 - term1 + R);
            stressYY = fx * (  - term2 - R);
            stressXY = fx * (  - term3 + S);
            stressYX = stressXY;
            stress = [ stressXX, stressXY; stressYX, stressYY ];
        end
        
        function volStressField = volStressField(mesh)
            % volStressField Analytical volume stress field
            
            nElements = mesh.numberOfElements;
            internalDataField = zeros(3,3,nElements);
            
            for iElement=1:nElements
                
                element = mesh.elements(iElement);
                internalDataField(:,:,iElement) = ...
                    demirdzic.PlateHole.stressElem(element);
            end
            internalField = tensorField(internalDataField);
            
            nBoundaries = mesh.numberOfBoundaries();
            bField = bTensorField(nBoundaries);
            
            for iBoundary = 1:nBoundaries
                
                boundary = mesh.boundaries(iBoundary);
                
                if strcmp(boundary.type,'empty') == true
                    
                    pField = patchTensorField(mesh, iBoundary);
                else
                    pDataField = zeros(3,3,boundary.numberOfBFaces);
                
                    startFace = boundary.startFace;
                    endFace = startFace + boundary.numberOfBFaces - 1;

                    for iFace=startFace:endFace

                        face = mesh.faces(iFace);
                        idx = iFace-startFace + 1;
                        pDataField(:,:,idx) = ...
                            demirdzic.PlateHole.stressElem(face);
                    end

                    data = tensorField(pDataField);
                    pField = ...
                        patchTensorField(mesh, iBoundary, data);
                end
                
                bField = bField.setPatch(iBoundary, pField);
            end
            
            volStressField = ...
                volTensorField(mesh, internalField, bField);
        end
        
        function plotContours(volStressField, nDivisions)
            % plotContours Plot fig. 9 [4.3; Demirdzic - 1994].
            
            option = struct(...
                'nDivXX',  1000:2000:29000,...
                'nDivYY', -9500:1000:4500,...
                'nDivXY', -8500:1000:1500....
                ); % Contour lines
            
            figHandle = figure(1);
            plotStress(...
                volStressField, nDivisions, [], [1 1], option.nDivXX);
            
            set...
            (...
                figHandle,...
                'NumberTitle', 'off', ...
                'Name', ...
                'sigma_xx  [fig. 9a, 4.3; Demirdzic - 1994]'...
            );

            figHandle = figure(2);
            plotStress(...
                volStressField, nDivisions, [], [2 2], option.nDivYY);
            
            set...
            (...
                figHandle,...
                'NumberTitle', 'off', ...
                'Name', ...
                'sigma_yy  [fig. 9b, 4.3; Demirdzic - 1994]'...
            );

            figHandle = figure(3);
            plotStress(...
                volStressField, nDivisions, [], [1 2], option.nDivXY);
            
            set...
            (...
                figHandle,...
                'NumberTitle', 'off', ...
                'Name', ...
                'sigma_xy  [fig. 9c, 4.3; Demirdzic - 1994]'...
            );
        end
    end
end

