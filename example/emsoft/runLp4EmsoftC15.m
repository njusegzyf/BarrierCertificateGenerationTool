function [lpVer, solveResVer, resNorms] = runLp4EmsoftC15()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftC15Problem();

x1 = vars(1);
x2 = vars(2);

degree = 2;
pLambdaDegree = 0;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = -45.9413122+45.5632*x1+72.6040*x2+26.2659*x1*x2-11.6836*x1^2-18.3213*x2^2;
% when add eps2 to theta, min rou : 1.7971e-05



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
