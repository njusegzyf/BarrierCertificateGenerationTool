function [lpVer, solveResVer, resNorms] = runLp4EmsoftDoub1()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftDoub1Problem();

x1 = vars(1);
x2 = vars(2);

degree = 2;
pLambdaDegree = 1;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = 120.5197041-104.9906*x1-37.8830*x2+41.0652*x1^2+2.1159*x2^2-25.5826*x1*x2;
% when add eps2 to theta, for pLambdaDegree = 1, min rou : 1.6332e-05
% when add eps2 to theta, for pLambdaDegree = 0, min rou : 1.7307e-05



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
