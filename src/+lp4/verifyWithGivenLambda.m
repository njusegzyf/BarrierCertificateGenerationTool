function [lpVer, solveResVer, resNorms] = verifyWithGivenLambda(...
    vars, f, eps, theta, psy, zeta, lambda, phyDegree)
%VERIFYWITHGIVENLLAMBDA 


import lp4.LinearProgram4Verification2
lp = LinearProgram4Verification2(vars);

lp.f = f;
lp.eps = eps;

% set the degree of phy and lambda
lp = lp.setDegreeAndInit(phyDegree);
lp.lambda = lambda;

lp = lp.setThetaConstraint(theta);
lp = lp.setPsyConstraint(psy);
lp = lp.setZetaConstraint(zeta);
lp = lp.generateEqsForConstraint1To3();

lp = lp.setDevVarsConstraint();

% solve the lp problem
[lpVer, solveResVer, resNorms] = lp.solve();

end
