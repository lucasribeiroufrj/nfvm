function checkBasicPossibleOperations(rawFields, fields)

mesh = fields{1}.mesh();

exps = ...
    {...
        @(lhs, rhs) sprintf('%s + %s', lhs, rhs),...
        @(lhs, rhs) sprintf('%s - %s', lhs, rhs),...
        @(lhs, rhs) sprintf('%s * %s', '0.87', rhs),...
        @(lhs, rhs) sprintf('%s * %s', lhs, '0.67'),...
        @(lhs, rhs) sprintf('%s / %s', lhs, '2.45')...
    };

for e = 1:length(exps)

    for i = 1:length(fields)

        field = fields{i};
        rawField = rawFields{i};   %#ok<NASGU>
        
        if strFind(class(field), 'Scalar')
            type = 'Scalar';
        elseif strFind(class(field), 'Vector')
            type = 'Vector';
        elseif strFind(class(field), 'Tensor')
            type = 'Tensor';
        end
        
        lmbd = exps{e};
        field = eval( lmbd('field','field') );
        rawField = eval( lmbd('rawField','rawField') );

        %% BEG - Check Internal rawField
        %
        isSurface = strFind(class(field), 'surface');
        
        if isSurface
            nElements = mesh.numberOfInteriorFaces;
        else
            nElements = mesh.numberOfElements;
        end
        
        len = min([length(size(rawField.data)) size(rawField.data) ]);
        switch len
            case 1
                data = rawField.data(1:nElements);
            case 2
                data = rawField.data(:,1:nElements);
            case 3
                data = rawField.data(:,:,1:nElements);
        end
        
        if ~all(all(field.internalField(1:nElements) == data))

            error('Test failed !');
        end
        %% END


        %% Check Boundary rawField
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

            len = min([length(size(rawField.data)) size(rawField.data) ]);
            switch len
                case 1
                    data = rawField.data(startIndex:endIndex);
                case 2
                    data = rawField.data(:,startIndex:endIndex);
                case 3
                    data = rawField.data(:,:,startIndex:endIndex);
            end

            patch = ...
                eval(['field.b',type,'Field_.patch',type,'Fields_{b}']);
            equality = data == patch(1:end);

            if ~all(all(equality))

                error('Test failed !');
            end
        end
    end
end

disp...
(...
    ['Basic possible operations. ',...
     '[PASSED]'...
    ]...
);

end
