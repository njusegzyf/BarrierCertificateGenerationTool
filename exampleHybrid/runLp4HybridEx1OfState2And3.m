function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4HybridEx1OfState2And3()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, stateNum, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards] = getLp4HybridEx1OfState2And3Problem();

% Set the degree
degree = 2;

res = [1, 1];
lambdas = [0, 0];

import lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe
lpVer = HybridLinearProgramVerificationWithGivenLambdaAndRe.create(vars, stateNum,...
    fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, degree,...
    lambdas, res);

[lpVer, solveResVer, resNorms] = lpVer.solve();
if solveResVer.hasSolution()
    return
end

warning('on')

echo off;
end
