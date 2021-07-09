function material = createMaterial(materialLawEnum, planeStress, varargin)
%createMaterial Summary of this function goes here
%   Detailed explanation goes here

materialLawName = char(materialLawEnum);

constructor = sprintf('materials.( ''%s'' )', materialLawName);

materialLawFactory = createObjectFactory(constructor);

material = materialLawFactory(planeStress, varargin{:});

end

