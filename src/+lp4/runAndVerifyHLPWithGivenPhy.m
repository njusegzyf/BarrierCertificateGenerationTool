function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyHLPWithGivenPhy(...
    lp, phyRange, pLambdaRange, pReRange)

lp = lp.setRangeAndWConstraint(phyRange, pLambdaRange, pReRange);

% solve the lp problem
[lp, solveRes] = lp.solve();

isVerified = false;

% verify the lp problem with the computed lambda
[lpVer, solveResVer, resNorms] = lp.verifyWithPhy(solveRes);
import lp4.isResNormsOk
if (solveRes.hasSolution() && solveResVer.hasSolution() && isResNormsOk(resNorms))
    isVerified = true;
    return;
end

end
