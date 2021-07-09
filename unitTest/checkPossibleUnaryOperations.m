function checkPossibleUnaryOperations(field, vField, types, exps) %#ok<INUSL>

mesh = vField.mesh();

for e = 1:length(exps)

    type = types{e};
    lmbd = exps{e};
    vTemp = eval( lmbd('vField') );
    temp = eval( lmbd('field') );

    %% BEG - Check Internal field
    %
    isSurface = strFind(class(vTemp), 'surface');
        
    if isSurface
        nElements = mesh.numberOfInteriorFaces;
    else
        nElements = mesh.numberOfElements;
    end

    if isa(temp,'scalarField')
        
            data = temp.data(1:nElements);
            
    elseif isa(temp, 'vectorField')
        
            data = temp.data(:,1:nElements);
            
    elseif isa(temp, 'tensorField')
        
            data = temp.data(:,:,1:nElements);
    end

    if ~all(all(vTemp.internalField(1:nElements) == data))

        error('Test failed !');
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

        if isa(temp,'scalarField')
        
            data = temp.data(startIndex:endIndex);
            
        elseif isa(temp, 'vectorField')

                data = temp.data(:,startIndex:endIndex);

        elseif isa(temp, 'tensorField')

                data = temp.data(:,:,startIndex:endIndex);
        end

        patch = eval(['vTemp.b',type,'Field_.patch',type,'Fields_{b}']);
        equality = data == patch(1:end);

        if ~all(all(equality))

            error('Test failed !');
        end
    end
    % END
end

disp...
(...
    ['Possible unary operations. ',...
     '[PASSED]'...
    ]...
);

end
