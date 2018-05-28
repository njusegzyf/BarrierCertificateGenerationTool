function [lpVer, solveResVer, resNorms] = runLp4EmsoftTest1()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftTest1Problem();

x = vars(1);
y = vars(2);

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

degree = 2;
pLambdaDegree = 1;

initPhy = x^2 + y^2 - 3;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
