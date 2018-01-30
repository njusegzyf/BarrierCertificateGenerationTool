function [lpVer, solveResVer, resNorms] = runLp4HybridEx1VerificationWithGivenLambdaAndRe()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx1Problem();
T = vars(1);
t = vars(2);

% Set the degree
degree = 2;

res = [1, 1, 1, 1];
lambdas = [0, 0, 0];

import lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe
lpVer = HybridLinearProgramVerificationWithGivenLambdaAndRe.create(vars, stateNum,...
    fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree,...
    lambdas, res);

[lpVer, solveResVer, resNorms] = lpVer.solve();
lp4.Lp4Config.displaySolveRes(solveResVer, resNorms);

% for lambda1 = [1, 0, -1]
%     for lambda2 = [1, 0, -1]
%         for lambda3 = [1, 0, -1]
%             
%             import lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe
%             lpVer = HybridLinearProgramVerificationWithGivenLambdaAndRe.create(vars, stateNum,...
%                 fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, degree,...
%                 [lambda1, lambda2, lambda3], res);
%             
%             [lpVer, solveResVer, resNorms] = lpVer.solve();
%             if solveResVer.hasSolution()
%                 return
%             end
%         end
%     end
% end

warning('on')
echo off;

end
