function [lpVer, solveResVer, resNorms] = verifyWithPhy(lp, solveRes)

resNorms = 0;
lpVer = 0;
solveResVer = 0;

if ~(solveRes.hasSolution())
    disp('No solution to verify.');
    return;
end

% verify

import lp4.LinearProgram4Verification3
lpVer = LinearProgram4Verification3(lp.indvars);

lpVer.f = lp.f;
lpVer.eps = lp.eps;

% Set the degree of lambda
lpVer = lpVer.setDegreeAndInit(lp.pLambdaDegree);

lpVer.phy = solveRes.getPhyExpression();

lpVer = lpVer.setThetaConstraint(lp.theta);
lpVer = lpVer.setPsyConstraint(lp.psy);
lpVer = lpVer.setZetaConstraint(lp.zeta);
lpVer = lpVer.generateEqsForConstraint1To3();

lpVer = lpVer.setDevVarsConstraint();

[lpVer, solveResVer, resNorms] = lpVer.solve();

if ~(solveResVer.hasSolution())
    disp('Verify failed.');
else
    disp('Verify succeed, norms ;');
    disp(resNorms);
end

end

