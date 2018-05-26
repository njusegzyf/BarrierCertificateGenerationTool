function [lpVer, solveResVer, resNorms] = runLp4Emsoft6d2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4Emsoft6d2Problem();

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);
x5 = vars(5);
x6 = vars(6);

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

degree = 2;
pLambdaDegree = 1;

initPhy = - 0.1182- 0.06094*x1- 0.00726*x2+ 0.3824*x3+ 0.3248*x4-0.1632*x5-0.1656*x6+ 0.03366*x1*x2- 0.006225*x1*x3- 0.004584*x1*x4- 0.009548*x1*x5+ 0.01588*x1*x6- 0.003326*x2*x3+ 0.005185*x2*x4+ 0.01407*x2*x5- 0.01693*x2*x6+ 0.04518*x3*x4  - 0.0167*x3*x5- 0.008248*x3*x6+ 0.004361*x4*x5+ 0.01167*x4*x6+ 0.0214*x5*x6  - 0.009488*x1^2- 0.009429*x2^2- 0.06402*x3^2- 0.08304*x4^2+0.005361*x5^2-0.007063*x6^2;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
