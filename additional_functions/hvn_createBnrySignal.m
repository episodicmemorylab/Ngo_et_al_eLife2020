% Function: hvn_createBnrySignal
% Created by H.-V.V. Ngo
%
% Creates a logical vector with '1'-bouts given by inBinary
%
% Usage: outBnry = hvn_createBnrySignal(inBnry, inLen)
%
% Input parameters
% inBnry    = [N x 2] array contraining N 'active'-bouts, with their on-
%             and offset depicted by the 1st and 2nd column, respectively
% inLen     = Length of resulting vector in data points
%%%

function outBnry = hvn_createBnrySignal(inBnry, inLen)

outBnry = zeros(1,inLen,'logical');

if ~isempty(inBnry)
    binaryIdx = cell2mat(reshape(arrayfun(@(x,y)  x:y, inBnry(:,1), inBnry(:,2),'UniformOutput',0),1,size(inBnry,1)));
    
    outBnry(binaryIdx) = 1;
end

end
