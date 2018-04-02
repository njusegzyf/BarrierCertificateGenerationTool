function [lpVer, solveResVer, resNorms] = runLp4Emsoft4d2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4Emsoft4d2Problem();

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);

degree = 2;
pLambdaDegree = 1;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = 5.782725478+34.8916*x1-72.8435*x2-3.6512*x3+70.1655*x4+5.1731*x1^2-4.3645*x2^2+2.3811*x3^2+4.6440*x2*x3-3.4586*x4^2-8.9548*x1*x2-2.9863*x1*x3+8.6539*x1*x4+8.6151*x2*x4-9.5623*x3*x4;
% when add eps2 to theta, pLambdaDegree = 1, min rou : 4.6644e-05



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
