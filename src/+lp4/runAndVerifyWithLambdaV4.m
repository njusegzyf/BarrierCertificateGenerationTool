function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithLambdaV4(...
    lp, phyRange, pLambdaRange, phyRangeInVerify)

lp.pPartitions = repmat(phyRange, 1024, 1);
lp.pLambdaPartitions = repmat(pLambdaRange, 1024, 1);
lp = lp.setWConstraint();

% solve the lp problem
[lp, solveRes] = lp.solve();

isVerified = false;

% verify the lp problem
[lpVer, solveResVer, resNorms] = lp.verify(solveRes, phyRangeInVerify);
import lp4.isResNormsOk
if (solveRes.hasSolution() && solveResVer.hasSolution() && isResNormsOk(resNorms))
    isVerified = true;
    return;
end

% if ~(solveResVer.hasSolution())
%     % split and run with bisected regions
% else
%     disp('Verify succeed, norms ;');
%     disp(resNorms);
%     return;
% end

end
