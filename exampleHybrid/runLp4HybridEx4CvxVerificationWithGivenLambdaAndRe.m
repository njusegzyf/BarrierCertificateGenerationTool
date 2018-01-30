function [lpVer, solveResVer, resNorms] = runLp4HybridEx4CvxVerificationWithGivenLambdaAndRe()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx4Problem();
x1 = vars(1);
x2 = vars(2);
x3 = vars(3);

% Set the degree
degree = 2;

res = [1, 1, 1, 1];
lambdas = [0, 0, 0];

import lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe
lpVer = HybridLinearProgramVerificationWithGivenLambdaAndRe.createWithRou(vars, stateNum,...
    fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree,...
    lambdas, res);

[lpVer, solveResVer, resNorms] = lpVer.solveWithCvx();

end
