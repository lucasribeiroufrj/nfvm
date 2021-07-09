function field = tryGetOrDefault(input, field, default)
% tryGetOrDefault(input, field, default) return input.field or default
%   depending if there is such field or not respectively.

    if isfield(input, field)
        
        field = input.(field);
    else
        field = default;
    end
end
