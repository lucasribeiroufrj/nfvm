function newVectorField = createVectorFieldFromInternalField(fileDirectory)
%createVectorFieldFromInternalField Summary of this function goes here
%   Detailed explanation goes here
%
% Dependency:
% - IODictionary
% - vectorField

dictionary = IODictionary(fileDirectory);
newVectorField = vectorField(dictionary('internalField'));

end

