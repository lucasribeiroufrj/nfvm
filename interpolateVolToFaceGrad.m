function surfaceGraPhi = interpolateVolToFaceGrad(...
    volGradPhi, phi)
%interpolateElementGradToFaceGrad Computes gradients at the face.
%   Interpolates gradients from volume field to surface field without 
% considering non-conjunctional mesh [Pag. 289 - Moukalled].

%% BEG - Set internal face gradient
mesh = volGradPhi.mesh;
wfField = mesh.wf();
wcField = mesh.wc();
dCFnormField = mesh.dCFnorm();
eCFfield = mesh.eCF();

iOwners = mesh.iOwners();
iNeighs = mesh.iNeighs();

wf = wfField.internalField();
wc = wcField.internalField();

gradC = tensorField(volGradPhi.internalField(iOwners));
gradF = tensorField(volGradPhi.internalField(iNeighs));
overBarGrad = wc*gradC + wf*gradF;

phiC = vectorField(phi.internalField(iOwners));
phiF = vectorField(phi.internalField(iNeighs));
dCFnorm = dCFnormField.internalField();
diffGrad = (phiF - phiC)/dCFnorm;

eCF = eCFfield.internalField();
K = diffGrad - overBarGrad*eCF;
M = zeros(3,3,mesh.numberOfInteriorFaces);
%M = tensorField(mesh.numberOfInteriorFaces);

L = scalarField(K.data(1,:))*eCF;
M(1,1,:) = L.data(1,:);
M(1,2,:) = L.data(2,:);
M(1,3,:) = L.data(3,:);

L = scalarField(K.data(2,:))*eCF;
M(2,1,:) = L.data(1,:);
M(2,2,:) = L.data(2,:);
M(2,3,:) = L.data(3,:);

L = scalarField(K.data(3,:))*eCF;
M(3,1,:) = L.data(1,:);
M(3,2,:) = L.data(2,:);
M(3,3,:) = L.data(3,:);

iFaceGraPhiField = tensorField(overBarGrad.data + M);
%% END

%% BEG - Set boundary face gradient
mesh = volGradPhi.mesh;

% "Mount" the surface field using internal and boundary fields.
surfaceGraPhi = ...
    surfaceTensorField(mesh, iFaceGraPhiField, volGradPhi.bTensorField());
%% END

end
