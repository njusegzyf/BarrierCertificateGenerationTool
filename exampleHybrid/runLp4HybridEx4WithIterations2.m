function [lpVer, solveResVer, resNorms] = runLp4HybridEx4WithIterations2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx4Problem();

% Set the degree
degree = 2;
pLambdaDegree = 2;
pReDegree = 2;
x1 = vars(1);
x2 = vars(2);
x3 = vars(2);

maxIterations = 50;

initRes = [1, 1];
initLambdas = [-0.2, -1];


import lp4.runAndVerifyHLPWithIterations2
[lpVer, solveResVer, resNorms] = runAndVerifyHLPWithIterations2(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initLambdas, initRes);

warning('on')
echo off;

end
