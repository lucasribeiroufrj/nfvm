function checkBoundaryIndexing(fields, nBoundaries)

types = {'Scalar', 'Vector', 'Tensor'};

kind = {'vol', 'surface'};

for k = 1:length(kind)
    
    kk = kind{k};
    
    for t = 1:length(types)
        
        tt = types{t};
        
        bField = eval(['fields{t}.b',tt,'Field_']); %#ok<NASGU>
        pField = eval(['bField.patch',tt,'Fields_']);
        
        for i = 1:nBoundaries
            
            switch t
                case 1
                    data = pField{i}.data(1:2);
                case 2
                    data = pField{i}.data(:,1:2);
                case 3
                    data = pField{i}.data(:,:,1:2);
            end
            
            equality = fields{t}.boundaryField(i).field(1:2) ...
                == data;

            if ~all(all(equality))

                error('Test failed !');
            end
        end
    end
end

disp...
(...
    ['Testing boundary indexing. ',...
     '[PASSED]'...
    ]...
);

end
