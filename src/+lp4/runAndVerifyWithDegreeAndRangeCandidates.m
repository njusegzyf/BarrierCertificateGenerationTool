function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithDegreeAndRangeCandidates(...
    vars, f, eps, theta, psy, zeta,...
    degrees, pLambdaDegrees,...
    phyRanges, pLambdaRanges, phyRangesInVerify)

rangeCandicatesCount = length(phyRanges);
if length(pLambdaRanges) ~= rangeCandicatesCount
    error('Inconsistent length of range candidates.');
end
if ~isempty(phyRangesInVerify) && length(phyRangesInVerify) ~= rangeCandicatesCount
    error('Inconsistent length of range candidates.');
end

isVerified = false;

for degree = degrees
    for pLambdaDegree = pLambdaDegrees
        import lp4.runAndVerifyWithRangeCandidatesV2
        [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithRangeCandidatesV2(...
            vars, f, eps, theta, psy, zeta,...
            degree, pLambdaDegree,...
            phyRanges, pLambdaRanges, phyRangesInVerify);
        
        if isVerified
            return;
        end
    end
end

end
