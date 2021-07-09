function checkPossibleBinaryOperations(rawFields, vFields, types, exps)

mesh = vFields{1}.mesh();

for e = 1:length(exps)

    v_lhs = vFields{1}; %#ok<NASGU>
    v_rhs = vFields{2}; %#ok<NASGU>
    f_lhs = rawFields{1};  %#ok<NASGU>
    f_rhs = rawFields{2};  %#ok<NASGU>
    type = types{e};

    lmbd = exps{e};
    vField = eval( lmbd('v_lhs','v_rhs') );
    rawField = eval( lmbd('f_lhs','f_rhs') );

    %% BEG - Check Internal field
    %
    isSurface = strFind(class(vField), 'surface');
        
    if isSurface
        nElements = mesh.numberOfInteriorFaces;
    else
        nElements = mesh.numberOfElements;
    end

    if isa(rawField,'scalarField')

        data = rawField.data(1:nElements);
            
    elseif isa(rawField, 'vectorField')

        data = rawField.data(:,1:nElements);
            
    elseif isa(rawField, 'tensorField')

        data = rawField.data(:,:,1:nElements);
    end

    if ~all(all(vField.internalField(1:nElements) == data))

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

        if isa(rawField,'scalarField')
        
            data = rawField.data(startIndex:endIndex);
            
        elseif isa(rawField, 'vectorField')

                data = rawField.data(:,startIndex:endIndex);

        elseif isa(rawField, 'tensorField')

                data = rawField.data(:,:,startIndex:endIndex);
        end

        patch = eval(['vField.b',type,'Field_.patch',type,'Fields_{b}']);
        equality = data == patch(1:end);

        if ~all(all(equality))

            error('Test failed !');
        end
    end
    % END
end

disp...
(...
    ['Basic possible operations. ',...
     '[PASSED]'...
    ]...
);

end
