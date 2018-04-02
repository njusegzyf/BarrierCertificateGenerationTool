function [lpVer, solveResVer, resNorms] = runLp4EmsoftC13_2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftC13_2Problem();

x1 = vars(1);
x2 = vars(2);

degree = 2;
pLambdaDegree = 0;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = 26.50266211-0.1484*x1-112.0672*x2+0.0508*x1*x2+45.2351*x1^2+5.4202*x2^2;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
