function [lpVer, solveResVer, resNorms] = runLp4EmsoftC6()
clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftC6Problem();

x1 = vars(1);
x2 = vars(2);

degree = 2;
pLambdaDegree = 0;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = 61.94063789+90.6204*x1+81.3145*x2-22.6008*x1^2-7.0696*x2^2+31.0851*x1*x2;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
