function checkEigensImplementation(field, vField, type)

mesh = vField.mesh();
[eigVetors, eigValues] = vField.eigen();
[eigVetorsF, eigValuesF] = field.eigen();

%% BEG - Check Internal field
%
isSurface = strFind(class(vField), 'surface');

if isSurface
    nElements = mesh.numberOfInteriorFaces;
else
    nElements = mesh.numberOfElements;
end

data1 = eigVetorsF.data(:,:,1:nElements);
data2 = eigValuesF.data(:,:,1:nElements);

if ~all(all(eigVetors.internalField(1:nElements) == data1)) | ...
   ~all(all(eigValues.internalField(1:nElements) == data2))

    error('Test failed!');
end
%% END


%% BEG - Check Boundary field
%
for b = 1:mesh.numberOfBoundaries

    boundary = mesh.boundaries(b);

    if isSurface
            startIndex = boundary.startFace;
    else
        startIndex = ....
            boundary.startFace - mesh.numberOfInteriorFaces...
            + mesh.numberOfElements;
    end

    nFaces = boundary.numberOfBFaces;
    endIndex = startIndex + nFaces - 1;

    data1 = eigVetorsF.data(:,:,startIndex:endIndex);
    data2 = eigValuesF.data(:,:,startIndex:endIndex);

    patch1 = eval(['eigVetors.b',type,'Field_.patch',type,'Fields_{b}']);
    patch2 = eval(['eigValues.b',type,'Field_.patch',type,'Fields_{b}']);
    equality1 = data1 == patch1(1:end);
    equality2 = data2 == patch2(1:end);

    if ~all(all(equality1)) | ~all(all(equality2))

        error('Test failed !');
    end
end
% END

disp...
(...
    ['Eigen operations. ',...
     '[PASSED]'...
    ]...
);

end
