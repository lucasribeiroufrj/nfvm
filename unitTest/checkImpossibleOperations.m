function checkImpossibleOperations(fields, vFields)

msg = 'Operation not allowed. ';

for i = 1:length(fields)
    
    field = fields{i};
    vField = vFields{i}; 
    
    checkIfFails(@() ( vField + field ), msg);
    checkIfFails(@() ( vField - field ), msg);
    checkIfFails(@() ( field + vField ), msg);
    checkIfFails(@() ( field - vField ), msg);
    checkIfFails(@() ( field * vField ), msg);
    checkIfFails(@() ( vField * field ), msg);
    checkIfFails(@() ( vField / field ), msg);
    checkIfFails(@() ( field / vField ), msg);
    
    if isempty(strfind(class(vField),'ScalarField'))
        checkIfFails(@() ( 1 / vField ), msg);
    end
end

end
