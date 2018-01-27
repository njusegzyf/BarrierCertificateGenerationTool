function [lpVer, solveResVer, resNorms] = runLp4BenkEx6AsLp()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx6Problem();
x = vars(1);
y = vars(2);



phyDegree = 2;
lambda = -1; % -6997937753528645/9007199254740992;



import lp4.verifyWithGivenLambda
[lpVer, solveResVer, resNorms] = verifyWithGivenLambda(vars, f, eps,  g_theta, g_psy, g_zeta, lambda, phyDegree);

import lp4.Lp4Config
Lp4Config.displaySolveRes(solveResVer, resNorms);

warning('on')

echo off;
end



% Verification succeeded.
% The parameter setting:phy degree: 2; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function phy is:
%    1.0e-03 *
% 
%     0.8259         0    0.9247         0   -0.0013         0
% 
% --------------------------------------------------------------
% The function phy is:
% (4264198615540793*y)/4611686018427387904 - (1520852061571313*x*y)/1180591620717411303424 + 238060677628973/288230376151711744
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0451
% 
% --------------------------------------------------------------
% Verify succeed, norms ;
%    1.0e-06 *
% 
%     0.0000    0.9559    0.0000
% 
% 
% ans = 
% 
%   LinearProgram4Verification2 -  Ù–‘:
% 
%           indvars: [1°¡2 sym]
%            degree: 2
%               phy: [1°¡1 sym]
%           phySize: 6
%            lambda: -0.7769
%               eps: [1.0000e-04 1.0000e-04]
%                 f: [2°¡1 sym]
%           decvars: [1°¡246 sym]
%             exprs: [1°¡4 Constraint]
%     phyPolynomial: [1°¡1 lp4util.SymbolicPolynomial]
%          solution: []
%          c1Length: 15
%          c2Length: 210
%          c3Length: 15
%       pPartitions: []
