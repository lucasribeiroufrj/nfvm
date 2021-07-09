function loadPaths()
    
paths = ...
    {...
        'fields',...
        'linearSolvers',...
        'solidModel',...
        'utilities',...
        'fvc',...
        'fvm',...
        'mesh',...
    };

root = pwd;
recursive = true;

addPath(paths, root, recursive);

paths = ...
    {...
        '',...
        'problems',...
        ['thirdParty' filesep 'MATAMG' filesep 'bin']...
    };
recursive = false;

addPath(paths, root, recursive);

    function addPathRecusive(path)

        addpath(genpath(path));
    end

    function addPathNotRecusive(path)

        addpath(path);
    end

    function fullPath = generateFullPath(root, path)
        
        fullPath = [root,'/', path];
    end

    function addPath(paths, root, recursive)
        
        if recursive
            
            lambda = @addPathRecusive;
        else
            lambda = @addPathNotRecusive;
        end
        
        for i = 1:length(paths)
            
            fullPath = generateFullPath(root, paths{i});
            lambda(fullPath);
        end
    end
end
