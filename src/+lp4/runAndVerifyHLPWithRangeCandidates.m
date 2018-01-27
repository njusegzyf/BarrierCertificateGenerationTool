function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyHLPWithRangeCandidates(...
    vars, stateNum, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards,...
    degree, pLambdaDegree, pReDegree,...
    phyRanges, pLambdaRanges, pReRanges)

rangeCandicatesCount = length(phyRanges);
if length(pLambdaRanges) ~= rangeCandicatesCount
    error('Inconsistent length of lambda range candidates.');
end
if length(pReRanges) ~= rangeCandicatesCount
    error('Inconsistent length of re range candidates.');
end

isVerified = false;

% Note: In this version, we reuse a `lp4.HybridLinearProgram` for different ranges.
import lp4.HybridLinearProgram
lp = HybridLinearProgram.createWithoutRanges(vars, stateNum, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards,...
    degree, pLambdaDegree, pReDegree);

for i = 1 : rangeCandicatesCount
    phyRange = phyRanges(i);
    pLambdaRange = pLambdaRanges(i);
    pReRange = pReRanges(i);
    
    import lp4.runAndVerifyHLPWithGivenPhy
    [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyHLPWithGivenPhy(...
        lp, phyRange, pLambdaRange, pReRange);
    
    % FIXME
    % import lp4.Lp4Config
    % if Lp4Config.IS_VERIFY_WITH_PHY
    % else
    % end
    
    if isVerified
        return;
    end
    
    % if a large range has no solution, we can skip small ranges
    if ~solveRes.hasSolution()
        return;
    end
end

end
