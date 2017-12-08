function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithDegreeAndRangeCandidates(...
    vars, f, eps, theta, psy, zeta,...
    degrees, pLambdaDegrees,...
    phyRanges, pLambdaRanges, phyRangesInVerify)

rangeCandicatesCount = length(phyRanges);
if length(pLambdaRanges) ~= rangeCandicatesCount
    error('Inconsistent length of range candidates.');
end
if phyRangesInVerify ~= 0 && length(phyRangesInVerify) ~= rangeCandicatesCount
    error('Inconsistent length of range candidates.');
end

isVerified = false;

for degree = degrees
    for pLambdaDegree = pLambdaDegrees
        import lp4.runAndVerifyWithRangeCandidates
        [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithRangeCandidates(...
            vars, f, eps, theta, psy, zeta,...
            degree, pLambdaDegree,...
            phyRanges, pLambdaRanges, phyRangesInVerify);
        
        if isVerified
            return;
        end
    end
end

end
