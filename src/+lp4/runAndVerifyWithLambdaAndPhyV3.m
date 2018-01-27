function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithLambdaAndPhyV3(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    phyRange, pLambdaRange, phyRangeInVerify)

import lp4.LinearProgram4_v3.createLp
lp = createLp(vars, f, eps, theta, psy, zeta, degree, pLambdaDegree, phyRange, pLambdaRange);

% solve the lp problem
[lp, solveRes] = lp.solve();

isVerified = false;

if ~solveRes.hasSolution()
	disp('Not find feasible solution to verify.');
    lpVer = 0;
    solveResVer = 0;
    resNorms = [];
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

% if ~(solveResVer.hasSolution())
%     % split and run with bisected regions
% else
%     disp('Verify succeed, norms ;');
%     disp(resNorms);
%     return;
% end

end
