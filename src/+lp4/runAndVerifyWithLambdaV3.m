function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithLambdaV3(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    phyRange, pLambdaRange, phyRangeInVerify)

import lp4.LinearProgram4_v3.createLp
lp = createLp(vars, f, eps, theta, psy, zeta, degree, pLambdaDegree, phyRange, pLambdaRange);

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
