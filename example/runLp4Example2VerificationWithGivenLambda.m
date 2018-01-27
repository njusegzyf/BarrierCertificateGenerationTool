function [lp, solveRes, resNorms] = runLp4Example2VerificationWithGivenLambda()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example2Problem();
x1 = vars(1);
x2 = vars(2);



phyDegree = 4;
lambda = -1.875;



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
% The parameter setting:phy degree: 4; eps1: 1e-05; eps2: 1e-05
% --------------------------------------------------------------
% The coefficients of function phy is:
%   1 至 13 列
% 
%     0.0032    0.0025    0.0019   -0.0002    0.0010   -0.0009   -0.0003    0.0004    0.0002   -0.0000    0.0001    0.0000    0.0001
% 
%   14 至 15 列
% 
%     0.0000   -0.0000
% 
% --------------------------------------------------------------
% The function phy is:
% (2673837054751623*x1^4)/36893488147419103232 + (2566938931570159*x1^3*x2)/295147905179352825856 - (5574431017609667*x1^3)/18446744073709551616 + (4871369375499747*x1^2*x2^2)/73786976294838206464 + (6796573574177367*x1^2*x2)/18446744073709551616 - (7943483086046341*x1^2)/36893488147419103232 + (8995819380003309*x1*x2^3)/2361183241434822606848 + (7911117701030957*x1*x2^2)/36893488147419103232 + (8783016330877579*x1*x2)/9223372036854775808 + (5820599625044703*x1)/2305843009213693952 - (6479611358767485*x2^4)/1180591620717411303424 - (3468068130345497*x2^3)/147573952589676412928 - (1975092145260681*x2^2)/2305843009213693952 + (4311105718440309*x2)/2305843009213693952 + 7327126933837585/2305843009213693952
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0425
% 
% --------------------------------------------------------------
% Verify succeed, norms ;
%    1.0e-17 *
% 
%     0.2868    0.1449    0.1338
% 
% 
% ans = 
% 
%   LinearProgram4Verification2 - 属性:
% 
%           indvars: [1×2 sym]
%            degree: 4
%               phy: [1×1 sym]
%           phySize: 15
%            lambda: -1.8750
%               eps: [1.0000e-05 1.0000e-05]
%                 f: [2×1 sym]
%           decvars: [1×365 sym]
%             exprs: [1×4 Constraint]
%     phyPolynomial: [1×1 lp4util.SymbolicPolynomial]
%          solution: []
%          c1Length: 70
%          c2Length: 210
%          c3Length: 70
%       pPartitions: []
      