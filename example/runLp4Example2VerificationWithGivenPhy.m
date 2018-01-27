function [lp, solveRes, resNorms] = runLp4Example2VerificationWithGivenPhy()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example2Problem();
x1 = vars(1);
x2 = vars(2);



lambdaDegree = 1;
phy = (2673837054751623*x1^4)/36893488147419103232 + (2566938931570159*x1^3*x2)/295147905179352825856 - (5574431017609667*x1^3)/18446744073709551616 + (4871369375499747*x1^2*x2^2)/73786976294838206464 + (6796573574177367*x1^2*x2)/18446744073709551616 - (7943483086046341*x1^2)/36893488147419103232 + (8995819380003309*x1*x2^3)/2361183241434822606848 + (7911117701030957*x1*x2^2)/36893488147419103232 + (8783016330877579*x1*x2)/9223372036854775808 + (5820599625044703*x1)/2305843009213693952 - (6479611358767485*x2^4)/1180591620717411303424 - (3468068130345497*x2^3)/147573952589676412928 - (1975092145260681*x2^2)/2305843009213693952 + (4311105718440309*x2)/2305843009213693952 + 7327126933837585/2305843009213693952;



import lp4.LinearProgram4Verification3
lp = LinearProgram4Verification3(vars);

lp.f = f;
lp.eps = eps;

% set the degree of lambda and phy
lp = lp.setDegreeAndInit(lambdaDegree);
lp.phy = phy;

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
%     0.2904
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
% 
% runLp4Example2VerificationWithGivenPhy
% constraint theta is processed: 2018-01-08 19:51:31
% constraint psy is processed: 2018-01-08 19:52:03
% constraint zeta is processed: 2018-01-08 19:52:19
% 
% Optimal solution found.
% 
% Verification succeeded.
% The parameter setting:lambda degree: 1; eps1: 1e-05; eps2: 1e-05
% --------------------------------------------------------------
% The parameter setting:
% lambdaDegree: 1; eps1: 1e-05; eps2: 1e-05
% --------------------------------------------------------------
% The coefficients of function lambda is:
%    -1.8750
%          0
%          0
% 
% --------------------------------------------------------------
% The function lambda is:
% -15/8
%  
% --------------------------------------------------------------
% The computation time is:
%     0.3692
% 
% --------------------------------------------------------------
% Verify succeed, norms ;
%    1.0e-16 *
% 
%     0.0330    0.1793    0.0129
% 
% 
% ans = 
% 
%   LinearProgram4Verification3 -  Ù–‘:
% 
%              indvars: [1°¡2 sym]
%         lambdaDegree: 1
%                  phy: [1°¡1 sym]
%               lambda: [1°¡1 sym]
%           lambdaSize: 3
%                  eps: [1.0000e-05 1.0000e-05]
%                    f: [2°¡1 sym]
%              decvars: [1°¡353 sym]
%                exprs: [1°¡4 Constraint]
%     lambdaPolynomial: [1°¡1 lp4util.SymbolicPolynomial]
%             solution: []
%             c1Length: 70
%             c2Length: 210
%             c3Length: 70
