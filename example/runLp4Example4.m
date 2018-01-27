function [lp, solveRes] = runLp4Example4()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example4Problem();



% Set the degree of phy and lambda
degree = 2;
pLambdaDegree = 0;

import lp4util.Partition
phyRange = Partition(-1, 1);
pLambdaRange = Partition(-1, 1);
phyRangeInVerify = Partition(-10, 10);



import lp4.runAndVerifyWithLambda
[lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithLambda(...
vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree,...
phyRange, pLambdaRange, phyRangeInVerify);

warning('on')

echo off;
end
