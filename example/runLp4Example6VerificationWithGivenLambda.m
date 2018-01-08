function [lpVer, solveResVer, resNorms] = runLp4Example6VerificationWithGivenLambda()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example6Problem();
x1 = vars(1);
x2 = vars(2);
x3 = vars(3);



phyDegree = 4;
lambda = -1;
% - (5144427543276353*x3)/1152921504606846976 - 2350404954671423/2251799813685248;



import lp4.verifyWithGivenLambda
[lpVer, solveResVer, resNorms] = verifyWithGivenLambda(vars, f, eps,  g_theta, g_psy, g_zeta, lambda, phyDegree);

import lp4.Lp4Config
Lp4Config.displaySolveRes(solveResVer, resNorms);

warning('on')

echo off;
end
