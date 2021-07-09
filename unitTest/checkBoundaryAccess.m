function field = ...
    checkBoundaryAccess(mesh, rawInput, fullTypeName, type)

    field = eval([fullTypeName,'(mesh, rawInput)']);
    for b = 1:mesh.numberOfBoundaries

        boundary = mesh.boundaries(b);
        
        if ~strFind(fullTypeName, 'surface')
        
            startIndex = ....
                boundary.startFace - mesh.numberOfInteriorFaces...
                + mesh.numberOfElements;
        else
            startIndex = boundary.startFace;
        end
        
        nFaces = boundary.numberOfBFaces;
        endIndex = startIndex + nFaces - 1;

        patch = eval(['field.b',type,'Field_.patch',type,'Fields_{b}']);
        equality = rawInput(startIndex:endIndex) == patch(1:end);

        if ~all(all(equality))

            error('Test failed !');
        end
    end

    disp...
    (...
        ['Constructing a ',fullTypeName,' from a raw input. ',...
         '[PASSED]'...
        ]...
    );
end
