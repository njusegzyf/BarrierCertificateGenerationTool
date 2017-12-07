function [lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithLambda(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    phyRange, pLambdaRange, phyRangeInVerify)

import lp4.LinearProgram4_v2
lp = LinearProgram4_v2(vars);

lp.f = f;
lp.eps = eps;

lp = lp.setDegreeAndInit(degree, pLambdaDegree);

lp = lp.setThetaConstraint(theta);
lp = lp.setPsyConstraint(psy);
lp = lp.setZetaConstraint(zeta);
lp = lp.generateEqsForConstraint1To3();

import lp4util.Partition
lp.pPartitions = repmat(phyRange, 1024, 1);
lp.pLambdaPartitions = repmat(pLambdaRange, 1024, 1);
lp = lp.setWConstraint();

lp = lp.setDevVarsConstraint();

lp = lp.setLinprogF();

% solve the lp problem
[lp, solveRes] = lp.solve();
% [lp, solveRes] = lp.solve1And3();

resNorms = 0;
lpVer = 0;
solveResVer = 0;

if ~(solveRes.hasSolution())
    disp('No solution to verify.');
    return;
end

% verify

import lp4.LinearProgram4Verification2
lpVer = LinearProgram4Verification2(vars);

lpVer.f = f;
lpVer.eps = eps;

% Set the degree of phy
lpVer = lpVer.setDegreeAndInit(degree);

lpVer.lambda = solveRes.getPLmabdaExpression();

lpVer = lpVer.setThetaConstraint(theta);
lpVer = lpVer.setPsyConstraint(psy);
lpVer = lpVer.setZetaConstraint(zeta);
lpVer = lpVer.generateEqsForConstraint1To3();

if isa(phyRangeInVerify, 'lp4util.Partition')
    lpVer.pPartitions = repmat(phyRangeInVerify, 1024, 1);
    lpVer.setPhyConstraint();
end

lpVer = lpVer.setDevVarsConstraint();

[lpVer, solveResVer, resNorms] = lpVer.solve();

if ~(solveResVer.hasSolution())
    disp('Verify failed.');
else
    disp('Verify succeed, norms ;');
    disp(resNorms);
end

% if ~(solveResVer.hasSolution())
%     % split and run with bisected regions
% else
%     disp('Verify succeed, norms ;');
%     disp(resNorms);
%     return;
% end

end
