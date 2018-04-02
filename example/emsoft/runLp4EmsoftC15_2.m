function [lpVer, solveResVer, resNorms] = runLp4EmsoftC15_2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftC15_2Problem();

x1 = vars(1);
x2 = vars(2);

degree = 2;
pLambdaDegree = 1;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = -44.80572631+44.6479*x1+71.3698*x2+25.8684*x1*x2-11.3483*x1^2-18.2959*x2^2;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
