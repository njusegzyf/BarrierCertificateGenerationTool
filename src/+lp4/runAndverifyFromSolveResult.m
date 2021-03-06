function  [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyFromSolveResult(...
    lp, solveRes, phyCoffDelta, lambdaCoffDelta, phyRangeInVerify)
%runAndVerifyFromSolveResult Reset the ranges of the `lp` based on the `solveRes`, and then rerun the lp and verify it.

if ~solveRes.hasSolution()
    error('');
end

phyCoff = solveRes.getPhyCoefficient();
lambdaCoff = solveRes.getPLmabdaCoefficient();

import lp4util.createPartitionsFromCoefficient
phyCoffPartitions = createPartitionsFromCoefficient(phyCoff, phyCoffDelta);
lambdaCoffPartitions = createPartitionsFromCoefficient(lambdaCoff, lambdaCoffDelta);

% reuse the lp
lp.pPartitions = phyCoffPartitions;
lp.pLambdaPartitions = lambdaCoffPartitions;
lp = lp.setWConstraint();

% solve the new lp problem
[lp, solveRes] = lp.solve();

isVerified = false;

% verify the new lp problem
[lpVer, solveResVer, resNorms] = lp.verify(solveRes, phyRangeInVerify);
import lp4.isResNormsOk
if (solveRes.hasSolution() && solveResVer.hasSolution() && isResNormsOk(resNorms))
    isVerified = true;
    return;
end

end

