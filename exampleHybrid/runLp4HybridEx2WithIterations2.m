function [lpVer, solveResVer, resNorms] = runLp4HybridEx2WithIterations2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx2Problem();
x1 = vars(1);
x2 = vars(2);

% Set the degree
degree = 1;
pLambdaDegree = 1;
pReDegree = 1;

maxIterations = 100;

initRes = [0, 0];
initLambdas = [0, 0];


import lp4.runAndVerifyHLPWithIterations2
[lpVer, solveResVer, resNorms] = runAndVerifyHLPWithIterations2(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initLambdas, initRes);

warning('on')
echo off;

end
