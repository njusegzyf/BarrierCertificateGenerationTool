function res = solveBlpDirectlyWithFmincon(vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree)

lp = lp4.LinearProgram4_v3(vars);

lp.f = f;
lp.eps = eps;

lp = lp.setDegreeAndInit(degree, pLambdaDegree);

lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

import lp4util.Partition

lp.pPartitions = repmat(Partition(-1, 1), [lp.lambdaEnd, 1]);
lp.pLambdaPartitions = repmat(Partition(-1, 1),[lp.lambdaEnd, 1]);
lp = lp.setWConstraint();

lp = lp.setDevVarsConstraint();

lp = lp.setLinprogF();



% solve the lp problem
res = lp.solveDirectlyWithFmincon();

end

