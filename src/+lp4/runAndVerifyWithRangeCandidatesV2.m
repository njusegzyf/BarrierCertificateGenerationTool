function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithRangeCandidatesV2(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    phyRanges, pLambdaRanges, phyRangesInVerify)

rangeCandicatesCount = length(phyRanges);
if length(pLambdaRanges) ~= rangeCandicatesCount
    error('Inconsistent length of range candidates.');
end
if ~isempty(phyRangesInVerify) && length(phyRangesInVerify) ~= rangeCandicatesCount
    error('Inconsistent length of range candidates.');
end

isVerified = false;

% Note: In this version, we reuse a `lp4.LinearProgram4_v3` for different ranges.
import lp4.LinearProgram4_v3
lp = LinearProgram4_v3.createLpWithoutRanges(vars, f, eps, theta, psy, zeta, degree, pLambdaDegree);

for i = 1 : rangeCandicatesCount
    phyRange = phyRanges(i);
    pLambdaRange = pLambdaRanges(i);
    if isempty(phyRangesInVerify)
        phyRangeInVerify = 0;
    else
        phyRangeInVerify = phyRangesInVerify(i);
    end
    
    import lp4.Lp4Config
    if Lp4Config.IS_VERIFY_WITH_PHY
        import lp4.runAndVerifyWithLambdaAndPhyV4
        [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithLambdaAndPhyV4(...
            lp, phyRange, pLambdaRange, phyRangeInVerify);
    else
        import lp4.runAndVerifyWithLambdaV4
        [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithLambdaV4(...
            lp, phyRange, pLambdaRange, phyRangeInVerify);
    end
    
    
    if isVerified
        return;
    end
    
    % if a large range has no solution, we can skip small ranges
    if ~solveRes.hasSolution()
        return;
    end
end

end
