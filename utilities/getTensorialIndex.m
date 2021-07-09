function interval = getTensorialIndex( index )
%getTensorialIndex(i) return the tensorial index of a block coupled
%   element.

startIdx = (index-1)*3 + 1;
endIdx = startIdx + 2;
interval = startIdx:endIdx;

end
