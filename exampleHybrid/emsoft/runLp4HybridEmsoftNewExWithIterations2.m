function [lpVer, solveResVer, resNorms] = runLp4HybridEmsoftNewExWithIterations2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEmsoftNewExProblem();
x = vars(1);

% Set the degree
degree = 1;
pLambdaDegree = 1;
pReDegree = 1;

maxIterations = 100;

% initRes = [6618173333825559/562949953421312 - (3610193783072245*x)/9007199254740992, (2038505519564809*x)/72057594037927936 - 5059550082599429/4503599627370496];
% initLambdas = [(1459638834325705*x)/332306998946228968225951765070086144 + 6606728061270839/281474976710656, 0];

% initRes = [1, 1];
% initLambdas = [1, 1];

% initRes = [0, 0];
% initLambdas = [0, 0];

initRes = [x, x];
initLambdas = [x, x];


[lpVer, solveResVer, resNorms] = lp4.runAndVerifyHLPWithIterations2(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initLambdas, initRes);

warning('on')
echo off;

end
