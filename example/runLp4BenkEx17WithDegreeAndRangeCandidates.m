function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4BenkEx17WithDegreeAndRangeCandidates()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx17Problem();



% Set the degree of phy and lambda
degrees = [1, 2, 3, 4];
pLambdaDegrees = [1, 2, 3];

ranges = [1, 0.5, 0.3, 0.15, 0.1];
import lp4util.createRangeCandidates
[phyRanges, pLambdaRanges, phyRangesInVerify] = createRangeCandidates(ranges, ranges, 0);



% run and verify
import lp4.runAndVerifyWithDegreeAndRangeCandidates
[lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithDegreeAndRangeCandidates(...
    vars, f, eps, g_theta, g_psy, g_zeta,...
    degrees, pLambdaDegrees,...
    phyRanges, pLambdaRanges, phyRangesInVerify);

warning('on')

echo off;
end



% --------------------------------------------------------------
% The parameter setting:
% degree: 1; lambda degree: 2; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function phy is:
%     0.4734   -0.6357   -0.3112
% 
% --------------------------------------------------------------
% The function phy is:
% 266486134759349/562949953421312 - (1401424382769031*y)/4503599627370496 - (59753*x)/94000
%  
% --------------------------------------------------------------
% The coefficients of lambda is:
%      0     0     0     0     0     0
% 
% --------------------------------------------------------------
% The function lambda is:
% 0
%  
% --------------------------------------------------------------
% The rou is:
%    -0.1578
% 
% --------------------------------------------------------------
% The computation time is:
%     0.0475
% 
% --------------------------------------------------------------
% constraint theta is processed: 2018-01-07 21:51:13
% constraint psy is processed: 2018-01-07 21:51:14
% constraint zeta is processed: 2018-01-07 21:51:14
% 
% No feasible solution found.
% 
% Linprog stopped because no point satisfies the constraints.
% 
% Verification failed.
% The parameter setting:phy degree: 1; eps1: 0.0001; eps2: 0.0001
% Verify feasible solution failed.
% --------------------------------------------------------------
% constraint theta is processed: 2018-01-07 21:51:18
% constraint psy is processed: 2018-01-07 21:51:19
% constraint zeta is processed: 2018-01-07 21:51:20
% 
% Optimal solution found.
% 
% Verification succeeded.
% The parameter setting:lambda degree: 2; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The parameter setting:
% lambdaDegree: 2; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function lambda is:
%     0.9633
%     0.8632
%     1.7046
%          0
%    -0.0000
%     0.0801
% 
% --------------------------------------------------------------
% The function lambda is:
% (3887295646963969*x)/4503599627370496 + (959614895919813*y)/562949953421312 - (2832394405306453*x*y)/81129638414606681695789005144064 + (1443341780062227*y^2)/18014398509481984 + 4338154958682301/4503599627370496
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0106
% 
% --------------------------------------------------------------
% Verify feasible solution succeed, norms ;
%    1.0e-15 *
% 
%          0    0.3941    0.1570
% 
% --------------------------------------------------------------
% 
% ans = 
% 
%   LinearProgram4_v3 -  Ù–‘:
% 
%                           indvars: [1°¡2 sym]
%                            degree: 1
%                               phy: [1°¡1 sym]
%                               eps: [1.0000e-04 1.0000e-04]
%                                 f: [2°¡1 sym]
%                           decvars: [1°¡73 sym]
%                             exprs: [1°¡5 Constraint]
%                     phyPolynomial: [1°¡1 lp4util.SymbolicPolynomial]
%                     pLambdaDegree: 2
%                 pLambdaPolynomial: [1°¡1 lp4util.SymbolicPolynomial]
%     pLambdaNormalizedSymbolicVars: []
%                     wSymbolicVars: [3°¡6 sym]
%                       wExpression: [1°¡1 sym]
%                            rouVar: [1°¡1 sym]
%                       pPartitions: [1024°¡1 lp4util.Partition]
%                 pLambdaPartitions: [1024°¡1 lp4util.Partition]
%                          linprogF: [1°¡73 double]
%                             theta: [1°¡2 sym]
%                               psy: [1°¡2 sym]
%                              zeta: [1°¡2 sym]
%                            rouInc: 0
