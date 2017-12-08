function [phyRangeCandidates, pLambdaRangeCandidates, phyRangesCandidatesInVerify] =...
    createRangeCandidates(phyRanges, pLambdaRanges, phyRangesInVerify)

rangesLength = length(phyRanges);
if length(pLambdaRanges) ~= rangesLength
    error('Inconsistent length of range candidates.');
end
if phyRangesInVerify ~= 0 && length(phyRangesInVerify) ~= rangesLength
    error('Inconsistent length of range candidates.');
end

import lp4util.Partition
verstring = char(version);
if verstring(1) == '9'
    % for new matlab version like 2017b, just use arrayfun
    phyRangeCandidates = arrayfun(@(x) Partition(-x, x), phyRanges);
    pLambdaRangeCandidates = arrayfun(@(x) Partition(-x, x), pLambdaRanges);
    if phyRangesInVerify == 0
        phyRangesCandidatesInVerify = 0;
    else
        phyRangesCandidatesInVerify = arrayfun(@(x) Partition(-x, x), phyRangesInVerify);
    end
else
    % for old matlab version like 2014b, 'arrayfun' cause errors so we directly construct arrays
    
    % pre allocated spaces
    phyRangeCandidates = repmat(Partition(0, 0), 1, rangesLength);
    pLambdaRangeCandidates = repmat(Partition(0, 0), 1, rangesLength);
    if phyRangesInVerify == 0
        phyRangesCandidatesInVerify = 0;
    else
        phyRangesCandidatesInVerify = repmat(Partition(0, 0), 1, rangesLength);
    end
    
    % fill range candidates
    for i = 1 : rangesLength
        phyRange = phyRanges(i);
        phyRangeCandidates(i) = Partition(-phyRange, phyRange);
        pLambdaRange = pLambdaRanges(i);
        pLambdaRangeCandidates(i) = Partition(-pLambdaRange, pLambdaRange);
        if (phyRangesInVerify ~= 0)
            phyRangeInVerify = phyRangesInVerify(i);
            phyRangesCandidatesInVerify(i) = Partition(-phyRangeInVerify, phyRangeInVerify);
        end
    end
end

end

