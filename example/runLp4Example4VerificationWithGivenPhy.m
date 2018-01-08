function [lp, solveRes, resNorms] = runLp4Example4VerificationWithGivenPhy()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example4Problem();
x1 = vars(1);
x2 = vars(2);



lambdaDegree = 1;
phy = (394901773787959*x1^2)/17592186044416 - (59123371936539*x1*x2)/8796093022208 + (2775580214094343*x1)/17592186044416 + (886329413328759*x2^2)/35184372088832 - (2151314439437839*x2)/35184372088832 - 2164372048877747/4398046511104;



import lp4.LinearProgram4Verification3
lp = LinearProgram4Verification3(vars);

lp.f = f;
lp.eps = eps;

% set the degree of lambda and phy
lp = lp.setDegreeAndInit(lambdaDegree);
lp.phy = phy;

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
