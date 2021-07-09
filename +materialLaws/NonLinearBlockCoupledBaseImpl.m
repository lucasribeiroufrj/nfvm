classdef NonLinearBlockCoupledBaseImpl ...
        < materialLaws.NonLinearBlockCoupledMaterialLaw ...
        & handle
    %NonLinearBlockCoupledBaseImpl Partial implementation for a NLBC
    %   material law.
    %   One can subclass this class to implement a NLBC material law
    %   without having to implementing everything needed.
    
    properties
        
    end
    
    methods
        %% BEG - materialLaws.NLBCMaterialLaw methods implementation
        %
        function Td = computeSurfaceTd(obj, F, k)
            % computeSurfaceTd return M . n (the . at 2).
            % The deformation gradient parameter F must be a surface
            % field. The parameter k is the row used from F.
            
            mesh = F.mesh();
            n = mesh.nf();
            t = F.rowAsVector(k);
            
            % T's internalField
            iTd = obj.computeTdTensorField(...
                F.internalField(),...
                n.internalField(),...
                t.internalField()...
             );
            
            nBoundaries = mesh.numberOfBoundaries;
            bTd = bTensorField(nBoundaries);

            for iBoundary = 1:nBoundaries

                % Get F's and t's patches.
                pF = F.boundaryField(iBoundary).field;
                pn = n.boundaryField(iBoundary).field;
                pt = t.boundaryField(iBoundary).field;

                % Use the tensor field from pC to create a patch field.
                field = obj.computeTdTensorField(pF,pn,pt);
                pTd = patchTensorField(mesh, iBoundary, field);

                % Add the calculated patch to R's boundary field.
                bTd = bTd.setPatch(iBoundary, pTd);
            end

            % "Mount" the final T field using the parts: 
            % internalfield + boundaryField
            Td = surfaceTensorField(mesh, iTd, bTd);
        end
        %% END - materialLaws.NLBCMaterialLaw methods implementation
    end
    
    methods (Abstract)
       
        % computeTdTensorField returns Td = FC : nt
        %                                    (2,3)
        % where C is elasticity tensor, n is the normal face and t is any
        % vector. The inputs must be tensor, vector and vector fields 
        % respectively.
        Td = computeTdTensorField(obj, F, n, t)
    end
end

