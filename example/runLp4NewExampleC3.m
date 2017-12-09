function [lp, solveRes, lpVer, solveResVer] = runLp4NewExampleC3()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC3Problem();



% Set the degree of phy and lambda
degree = 6;
pLambdaDegree = 2;

import lp4util.Partition
phyRange = Partition(-1, 1);
pLambdaRange = Partition(-1, 1);
phyRangeInVerify = 0; % Partition(-10, 10);



import lp4.runAndVerifyWithLambdaV3
[lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithLambdaV3(...
vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree,...
phyRange, pLambdaRange, phyRangeInVerify);

warning('on')

echo off;
end
