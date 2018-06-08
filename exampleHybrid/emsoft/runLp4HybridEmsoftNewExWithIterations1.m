function [lpVer, solveResVer, resNorms] = runLp4HybridEmsoftNewExWithIterations1()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEmsoftNewExProblem();
x = vars(1);


% Set the degree
degree = 2;
pLambdaDegree = 2;
pReDegree = 2;

maxIterations = 20;

initPhys = [591232478353217/147573952589676412928 - (2314583338575937*x)/4722366482869645213696,...
            (8198255338659429*x)/1180591620717411303424 - 570617907208669/2305843009213693952];
initPhys = initPhys * 1e5;



[lpVer, solveResVer, resNorms] = lp4.runAndVerifyHLPWithIterations1(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initPhys);

warning('on')
echo off;

end
