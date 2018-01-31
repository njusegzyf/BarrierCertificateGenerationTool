function [lpVer, solveResVer, resNorms] = runLp4Example9WithIteration2()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, theta, psy, zeta] = getLp4Example9Problem();

degree = 2;
pLambdaDegree = 1;

maxIterations = 100;

initLambda = 0;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations2(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initLambda);

warning('on')

echo off;
end
