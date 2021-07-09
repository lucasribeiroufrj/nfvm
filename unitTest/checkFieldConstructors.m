function checkFieldConstructors(mesh) %#ok<INUSD>

types = {'Scalar', 'Vector', 'Tensor'};
kind = {'vol','surface'};

for k = 1:length(kind)
    
    kk = kind{k};
    
    for t = 1:length(types)
        
        tt = types{t};
        
        handle1 = eval(['@() ',kk,tt,'Field()']);
        handle2 = eval(['@() ',kk,tt,'Field(mesh)']);
        handle3 = eval(['@() ',kk,tt,'Field(mesh,mesh.numberOfElements)']);
        
        checkIfFails(handle1,...
            [kk,tt,'Field:Default construction should be impossible. '])

        checkIfPasses(handle2,...
            [kk,tt,'Field:Construction from mesh must be possible. '])

        checkIfFails(handle3,...
            [kk,tt,'Field:Construction from int should be impossible. '])
    end
end

end
