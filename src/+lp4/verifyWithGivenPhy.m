function [lpVer, solveResVer, resNorms] = verifyWithGivenPhy(...
    vars, f, eps, theta, psy, zeta, phy, lambdaDegree)
%VERIFYWITHGIVENPHY 

import lp4.LinearProgram4Verification3
lp = LinearProgram4Verification3(vars);

lp.f = f;
lp.eps = eps;

% set the degree of lambda and phy
lp = lp.setDegreeAndInit(lambdaDegree);
lp.phy = phy;

lp = lp.setThetaConstraint(theta);
lp = lp.setPsyConstraint(psy);
lp = lp.setZetaConstraint(zeta);
lp = lp.generateEqsForConstraint1To3();

lp = lp.setDevVarsConstraint();

% solve the lp problem
[lpVer, solveResVer, resNorms] = lp.solve();

end
