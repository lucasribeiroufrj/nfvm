classdef SegregatedBaseImpl ...
        < materialLaws.SegregatedMaterialLaw ...
        & handle
    %SegregatedBaseImpl Partial implementation for a Segregated material
    %   law.
    %   One can subclass this class to implement a Segregated material law
    %   without having to implementing everything needed.
    
    properties
       
        % Material properties
        material_;
        
        % Implicit stiffness
        surfaceImpK_;
        volImpK_;
    end
    
    methods
        
        function obj = SegregatedBaseImpl(mesh, material)
            % SegregatedBaseImpl Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.material_ = material;
            
            %% BEG - Set Impk
            
            data = ones(1,mesh.numberOfFaces) ...
                *(2*obj.material_.mu() + obj.material_.lambda());
            sImpK = scalarField(data);
            N = mesh.numberOfElements + mesh.numberOfBElements;
            vImpK = sImpK.subField(1:N);

            obj.surfaceImpK_ = surfaceScalarField(mesh, sImpK);
            obj.volImpK_ = volScalarField(mesh, vImpK);
            %% END - Set Impk
        end
        
        %% BEG - Data member access
        %
        function material = material(obj)
            
            material = obj.material_;
        end

        % Beg - materialLaws.MaterialLaw interface implementation
        function density = density(obj)
            
            density = obj.material().density();
        end
        
        function planeStress = planeStress(obj)
            
            planeStress = obj.material().planeStress();
        end
        % End - materialLaws.MaterialLaw interface implementation
        
        function surfaceImpK = surfaceImpK(obj)
            
            surfaceImpK = obj.surfaceImpK_;
        end
        
        function volImpK = volImpK(obj)
            
            volImpK = obj.volImpK_;
        end
        
        function obj = setSurfaceImpK(obj, surfaceImpK)
            
            obj.surfaceImpK_ = surfaceImpK;
        end
        
        function obj = setVolImpK(obj, volImpK)
            
            obj.volImpK_ = volImpK;
        end
        %% END - Data member acces
        
        
        %% BEG - SegregatedMaterialLaw interface implementation
        %
        function Piola = computePiolaField(obj, gradU)
            % computePiolaField the first-piola stress tensor for a given
            %   displacement gradient field. The gradUField parameter can
            %   be a volume or a surface fields.
            
            I = gradU.identity();
            F = I + gradU;

            Sigma = obj.computeSigmaField(gradU);
            Piola = F*Sigma;
        end
        %% END - SegregatedMaterialLaw interface implementation
        
        
        %% BEG - Other member functions
        %
        function Sigma = computeSigmaField(obj, gradU)
            % computeSigmaField Computes the second-piola stress 
            %   tensor for a given displacement gradient field.
            %   gradU : volTensorField, surfaceTensorField or 
            %   patchTensorField.
            
            if isa(gradU, 'patchTensorField')
                
                Sigma = obj.computeSigmaPatchField(gradU);
                return;
            end
            
            % Sigma's internal field.
            iSigma = obj.computeSigmaTensorField(gradU.internalField());
            
            nBoundaries = gradU.mesh().numberOfBoundaries;
            bSigmaField = bTensorField(nBoundaries);
            mesh = gradU.mesh();
            
            for iBoundary = 1:nBoundaries
                
                % Get gradU's patch.
                pGradU = gradU.boundaryField(iBoundary).field;
                
                % Use the patch to calculate a tensor field.
                SigmaField = obj.computeSigmaTensorField(pGradU);
                
                % Use the tensor field to create a patch field.
                pSigma = patchTensorField(mesh, iBoundary, SigmaField);
                
                % Add the calculated patch to Sigma's boundary field.
                bSigmaField = bSigmaField.setPatch(iBoundary, pSigma);
            end
            
            % "Mount" the final Sigma field using the parts: 
            % internalfield + boundaryField
            if isa(gradU, 'surfaceTensorField')
                
                Sigma = surfaceTensorField(mesh, iSigma, bSigmaField);
            else
                Sigma = volTensorField(mesh, iSigma, bSigmaField);
            end
        end
        
        function Sigma = computeSigmaPatchField(obj, gradU)
            % computeSigmaPatchField Computes the second-piola stress.
            %   gradU must be a patchField.
            
            SigmaField = obj.computeSigmaTensorField(gradU.field());
            
            Sigma = ...
                patchTensorField(gradU.mesh(), gradU.index(), SigmaField);
        end
        %% END - Other member functions
    end
    
    methods (Abstract)
        
        % Computes the second piola stress given the right-Cauchy-Green
        %   strain tansor (as a tensorField).
        %
        %   computeSigmaTensorField(gradU) returns a tensorField. The
        %       argument gradU must be tensorField too.
        Sigma = computeSigmaTensorField(obj, gradU)
    end
end

