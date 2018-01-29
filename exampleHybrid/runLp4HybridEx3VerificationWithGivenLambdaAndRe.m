function [lpVer, solveResVer, resNorms] = runLp4HybridEx3VerificationWithGivenLambdaAndRe()

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

res = [1, 1, 1, 1];
lambdas = [0, 0, 0];

import lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe
lpVer = HybridLinearProgramVerificationWithGivenLambdaAndRe.create(vars, stateNum,...
    fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree,...
    lambdas, res);

[lpVer, solveResVer, resNorms] = lpVer.solve();
if solveResVer.hasSolution() && lp4.isResNormsOk(resNorms)
    return;
end

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



% res = [1, 1, 1, 1];
% lambdas = [0, 0, 0];
% 
%The parameter setting:
% degree: 1; eps1: 0.01; eps2: 0.01
% --------------------------------------------------------------
% The coefficients of function phy1 is:
%   -0.326666666666667   0.062222222222222   0.073333333333333
% 
% The function phy1 is:
% (14*x1)/225 + (11*x2)/150 - 49/150
%  
% The coefficients of function phy2 is:
%    0.370000000000000  -0.050000000000000   0.050000000000000
% 
% The function phy2 is:
% x2/20 - x1/20 + 37/100
