function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4BenkEx6AsLp()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx6Problem();
x = vars(1);
y = vars(2);



phyDegree = 2;
lambda = -6997937753528645/9007199254740992;



import lp4.verifyWithGivenLambda
[lpVer, solveResVer, resNorms] = verifyWithGivenLambda(vars, f, eps,  g_theta, g_psy, g_zeta, lambda, phyDegree);

import lp4.Lp4Config
Lp4Config.displaySolveRes(solveResVer, resNorms);

warning('on')

echo off;
end
