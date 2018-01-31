function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithLambdaAndPhyV4(...
    lp, phyRange, pLambdaRange, phyRangeInVerify)

% Compared to `runAndVerifyWithLambdaAndPhyV3`, this version use(reuse) the given LP for verification.

lp.pPartitions = repmat(phyRange, 1024, 1);
lp.pLambdaPartitions = repmat(pLambdaRange, 1024, 1);
lp = lp.setWConstraint();

% solve the lp problem
[lp, solveRes] = lp.solve();

isVerified = false;

if ~solveRes.hasSolution()
    disp('Not find feasible solution to verify.');
    [lpVer, solveResVer, resNorms] = lp4.Lp4Config.createAbsentVerificationResult();
    return;
end

% verify the lp problem with the computed lambda
[lpVer, solveResVer, resNorms] = lp.verifyWithLambda(solveRes, phyRangeInVerify);
import lp4.isResNormsOk
if solveRes.hasSolution() && solveResVer.hasSolution() && isResNormsOk(resNorms)
    isVerified = true;
    return;
end

% verify the lp problem with the computed phy
[lpVer, solveResVer, resNorms] = lp.verifyWithPhy(solveRes);
import lp4.isResNormsOk
if solveRes.hasSolution() && solveResVer.hasSolution() && isResNormsOk(resNorms)
    isVerified = true;
    return;
end

end
