function [lpVer, solveResVer, resNorms] = runLp4Emsoft6d1()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4Emsoft6d1Problem();

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);
x5 = vars(5);
x6 = vars(6);

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

degree = 2;
pLambdaDegree = 1;

initPhy = - 0.1725*x1^2 - 0.6077*x1*x2 - 0.5707*x1*x3 - 0.03331*x1*x4 - 0.47*x1*x5 - 0.03331*x1*x6 + 0.3315*x1 + 0.209*x2^2 + 0.2005*x2*x3 - 0.03296*x2*x4 - 0.03367*x2*x5 - 0.03296*x2*x6 + 0.5343*x2 + 2.997*x3^2 - 0.2995*x3*x4 - 0.3037*x3*x5 - 0.2995*x3*x6 + 5.1*x3 + 1.427*x4^2 - 0.3047*x4*x5 - 0.3164*x4*x6 + 1.916*x4 + 1.989*x5^2 - 0.3047*x5*x6 + 2.912*x5 + 1.427*x6^2 + 1.916*x6 + 6.718;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
