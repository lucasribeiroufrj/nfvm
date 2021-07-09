function factory = createObjectFactory(className)
% createObjectFactory(className) return an anonymous constructor function
%   for type named className.
%   Ex.: 
%   >> fvMeshFactory = createObjectFactory('fvMesh');
%   >> a_timeMesh = ...;
%   >> a_spaceMesh = ...;
%   >> obj.fvMesh = fvMeshFactory(a_timeMesh, a_spaceMesh);
    
factory = eval(['@(varargin)', className,'(varargin{:});']);

end
