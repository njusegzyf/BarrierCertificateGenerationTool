function [lpVer, solveResVer, resNorms] = runLp4EmsoftC14()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftC14Problem();
x1 = vars(1);
x2 = vars(2);
x3 = vars(3);

degree = 2;
pLambdaDegree = 1;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

initPhy = 17.9226844-5.4569*x1-0.7915*x2-19.4404*x3+50.2870*x1^2-30.3634*x1*x2+4.8555*x1*x3+12.9877*x2^2-19.7978*x2*x3+0.7835*x3^2;
% when add eps2 to theta, min rou : 1.4286e-05



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
