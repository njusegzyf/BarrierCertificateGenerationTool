function [lpVer, solveResVer, resNorms] = runLp4EmsoftC6_2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftC6_2Problem();

x1 = vars(1);
x2 = vars(2);

degree = 2;
pLambdaDegree = 0;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = 69.46478566+94.0516*x1+83.3774*x2-20.3721*x1^2-6.6288*x2^2+32.0463*x1*x2;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
