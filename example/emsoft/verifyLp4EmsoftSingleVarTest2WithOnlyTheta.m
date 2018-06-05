function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftSingleVarTest2WithOnlyTheta()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, ~, ~] = getLp4EmsoftSingleVarTest2Problem();

x = vars(1);

initPhy = (2*x-1)^2;

pLambdaDegree = 0;



% run and verify
lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, [], [], pLambdaDegree, initPhy);
[lpVer, solveResVer, resNorms] = lpVer.solve();

lp4.Lp4Config.printVerifyWithOnlyPsyResult(solveResVer, resNorms);

warning('on')

echo off;
end
