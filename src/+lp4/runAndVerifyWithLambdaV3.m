function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithLambdaV3(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    phyRange, pLambdaRange, phyRangeInVerify)

% Note: Matlab do not support `import lp4.LinearProgram4_v3.createLp`
% (a static method) and then use it directly.
import lp4.LinearProgram4_v3
lp = LinearProgram4_v3.createLp(vars, f, eps, theta, psy, zeta, degree, pLambdaDegree, phyRange, pLambdaRange);

% solve the lp problem
[lp, solveRes] = lp.solve();

isVerified = false;

% verify the lp problem with the computed lambda
[lpVer, solveResVer, resNorms] = lp.verifyWithLambda(solveRes, phyRangeInVerify);
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
