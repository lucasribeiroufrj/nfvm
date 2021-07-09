function bool = strFind( str1, str2 )
%strFind Checks if one string contains the other
%   Detailed explanation goes here

if ~isempty(strfind(str1,str2))
    
    bool = true;
else
    bool = false; 
end

end

