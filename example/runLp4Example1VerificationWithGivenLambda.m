function [lp, solveRes, resNorms]  = runLp4Example1VerificationWithGivenLambda()

% Example 1 

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [-(11/2)*x2+x2^2;
    6*x1-x1^2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 2*x1-8;
theta2 = x2;
theta3 = x2-1;

g_theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = (x1-1)/4;
psy2 = (x2-1)/4;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1-1;
zeta2 = x2-2;

g_zeta = [zeta1, zeta2];



phyDegree = 3;
lambda = - (1125925798617485*x1)/2251799813685248 - (1125925798617485*x2)/2251799813685248 - 1125925798617485/2251799813685248;



import lp4.LinearProgram4Verification2
lp = LinearProgram4Verification2(vars);

lp.f = f;
lp.eps = eps;

% set the degree of phy and lambda
lp = lp.setDegreeAndInit(phyDegree);
lp.lambda = lambda;

lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

lp = lp.setDevVarsConstraint();

% solve the lp problem
[lp, solveRes, resNorms] = lp.solve();

if ~(solveRes.hasSolution())
    disp('Verify failed.');
else
    disp('Verify succeed, norms ;');
    disp(resNorms);
end

warning('on')

echo off;
end



% Verification succeeded.
% The parameter setting:phy degree: 3; eps1: 0.001; eps2: 0.001
% --------------------------------------------------------------
% The coefficients of function phy is:
%    -0.0411   -0.0081    0.0074    0.0050    0.0007    0.0025   -0.0002   -0.0002    0.0002   -0.0004
% 
% --------------------------------------------------------------
% The function phy is:
% - (2263286168729447*x1^3)/9223372036854775808 - (8341502826980283*x1^2*x2)/36893488147419103232 + (5739950446431231*x1^2)/1152921504606846976 + (7880846712646413*x1*x2^2)/36893488147419103232 + (6463301049956541*x1*x2)/9223372036854775808 - (4649760513280169*x1)/576460752303423488 - (7880665484799027*x2^3)/18446744073709551616 + (5749322513946143*x2^2)/2305843009213693952 + (8475741288162059*x2)/1152921504606846976 - 5925416967917561/144115188075855872
%  
% --------------------------------------------------------------
% The computation time is:
%     0.2494
% 
% --------------------------------------------------------------
% Verify succeed, norms ;
%    1.0e-16 *
% 
%     0.4420    0.9089    0.0664
% 
% 
% ans = 
% 
%   LinearProgram4Verification2 -  Ù–‘:
% 
%           indvars: [1°¡2 sym]
%            degree: 3
%               phy: [1°¡1 sym]
%           phySize: 10
%            lambda: [1°¡1 sym]
%               eps: [1.0000e-03 1.0000e-03]
%                 f: [2°¡1 sym]
%           decvars: [1°¡199 sym]
%             exprs: [1°¡4 Constraint]
%     phyPolynomial: [1°¡1 lp4util.SymbolicPolynomial]
%          solution: []
%          c1Length: 84
%          c2Length: 70
%          c3Length: 35
%       pPartitions: []
