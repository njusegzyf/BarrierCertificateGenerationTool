function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithRangeCandidates2(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    phyRanges, pLambdaRanges, phyRangesInVerify)

rangeCandicatesCount = length(phyRanges);
if length(pLambdaRanges) ~= rangeCandicatesCount
    error('Inconsistent length of range candidates.');
end
if phyRangesInVerify ~= 0 && length(phyRangesInVerify) ~= rangeCandicatesCount
    error('Inconsistent length of range candidates.');
end

isVerified = false;

for i = 1 : rangeCandicatesCount
    phyRange = phyRanges(i);
    pLambdaRange = pLambdaRanges(i);
    if phyRangesInVerify == 0
        phyRangeInVerify = 0;
    else
        phyRangeInVerify = phyRangesInVerify(i);
    end
    
    import lp4.runAndVerifyWithLambdaAndPhyV3
    [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithLambdaAndPhyV3(...
        vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
        phyRange, pLambdaRange, phyRangeInVerify);
    
     if isVerified
         return;
     end
end

end
