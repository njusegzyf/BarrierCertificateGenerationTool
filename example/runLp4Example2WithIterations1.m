function [lpVer, solveResVer, resNorms] = runLp4Example2WithIterations1()

% Ok in ROU_THRESHOLD = 1e-6

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4Example2Problem();
x1 = vars(1);
x2 = vars(2);

degree = 2;
pLambdaDegree = 1;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = (66853*x1*x2)/8796093022208 + (2117713*x1^2)/8796093022208 - (1394370808931*x2^2)/17592186044416;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
