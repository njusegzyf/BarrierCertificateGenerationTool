function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyHLPWithDegreeAndRangeCandidates(...
    vars, stateNum, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards,...
    degrees, pLambdaDegrees, pReDegrees,...
    phyRanges, pLambdaRanges, pReRanges)

rangeCandicatesCount = length(phyRanges);
if length(pLambdaRanges) ~= rangeCandicatesCount
    error('Inconsistent length of lambda range candidates.');
end
if length(pReRanges) ~= rangeCandicatesCount
    error('Inconsistent length of re range candidates.');
end

isVerified = false;

for degree = degrees
    for pLambdaDegree = pLambdaDegrees
        for pReDegree = pReDegrees
            import lp4.runAndVerifyHLPWithRangeCandidates
            [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyHLPWithRangeCandidates(...
                vars, stateNum, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards,...
                degree, pLambdaDegree, pReDegree,...
                phyRanges, pLambdaRanges, pReRanges);
            
            if isVerified
                return;
            end
        end
    end
end

end
