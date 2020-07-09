% Function: hvn_extrctBnryBouts
% Created by H.-V.V. Ngo
%
% Extracts on- and offset of 'active'-bouts within a logical vector
%
% Usage: outBnry = hvn_extrctBnryBouts(inBnry)
%
% Input parameters
% inBnry    = 1-D logical vector with active bouts (periods with value 1)
%
% Output
% outBnry   = [N x 2] vector of N bouts with on- and offset represented by
%             1st and 2nd column
%%%

function outBnry = hvn_extrctBnryBouts(inBnry)

if size(inBnry,1) ~= 1
    inBnry = inBnry';
end

bnryOn      = find(diff([0 inBnry]) ==  1);
bnryOff     = find(diff([inBnry 0]) == -1);
numBouts    = min(numel(bnryOn),numel(bnryOff));

outBnry = [bnryOn(1:numBouts)', bnryOff(1:numBouts)'];

end
