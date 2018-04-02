function [lpVer, solveResVer, resNorms] = runLp4EmsoftC14_2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftC14_2Problem();
x1 = vars(1);
x2 = vars(2);
x3 = vars(3);

degree = 2;
pLambdaDegree = 0;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = 165.1787867-115.3431*x1+16.1618*x2-210.1866*x3+51.0304*x1^2-23.1069*x1*x2+4.9340*x1*x3+9.9151*x2^2+31.0738*x2*x3-11.6747*x3^2;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
