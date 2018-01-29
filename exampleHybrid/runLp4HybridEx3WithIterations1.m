function [lpVer, solveResVer, resNorms] = runLp4HybridEx3WithIterations1()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx3Problem();
x1 = vars(1);
x2 = vars(2);

% Set the degree
degree = 1;
pLambdaDegree = 1;
pReDegree = 1;

maxIterations = 10;

initPhys = [(14*x1)/225 + (11*x2)/150 - 50/150,...
             x2/20 - x1/20 + 50/100];



import lp4.runAndVerifyHLPWithIterations1
[lpVer, solveResVer, resNorms] = runAndVerifyHLPWithIterations1(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initPhys);

warning('on')
echo off;

end
