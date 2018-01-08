function [lp, solveRes, resNorms] = runLp4Example2VerificationWithGivenLambda()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example2Problem();
x1 = vars(1);
x2 = vars(2);



phyDegree = 4;
lambda = -1.875;



import lp4.LinearProgram4Verification2
lp = LinearProgram4Verification2(vars);

lp.f = f;
lp.eps = eps;

% set the degree of phy and lambda
lp = lp.setDegreeAndInit(phyDegree);
lp.lambda = lambda;

lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

lp = lp.setDevVarsConstraint();

% solve the lp problem
[lp, solveRes, resNorms] = lp.solve();

if ~(solveRes.hasSolution())
    disp('Verify failed.');
else
    disp('Verify succeed, norms ;');
    disp(resNorms);
end

warning('on')

echo off;
end
