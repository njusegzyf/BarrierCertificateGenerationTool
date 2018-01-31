function [lpVer, solveResVer, resNorms] = runLp4Example3WithIterations1()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4Example3Problem();
x = vars(1);
y = vars(2);

degree = 2;
pLambdaDegree = 1;

maxIterations = 10;

initPhy = (66853*x*y)/8796093022208 + (2117713*x^2)/8796093022208 - (1394370808931*y^2)/17592186044416;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end


