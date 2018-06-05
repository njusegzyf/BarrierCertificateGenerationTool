function [lp, solveRes] = runLp4Example7()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example7Problem();



% Set the degree of phy and lambda
degree = 2;
pLambdaDegree = 0;

import lp4util.Partition
phyRange = Partition(-0.3, 0.3);
pLambdaRange = Partition(-0.3, 0.3);



import lp4.LinearProgram4_v3
lp = LinearProgram4_v3.createLp(vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree, phyRange, pLambdaRange);

% solve the lp problem
[lp, solveRes] = lp.solve();

echo off;
end
