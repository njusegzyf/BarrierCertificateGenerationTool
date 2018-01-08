function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4Example3WithDegreeAndRangeCandidates()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example3Problem();
                              

% set the degree of phy and lambda
degrees = [1,2];
pLambdaDegrees = [0,1,2,3];

% set the ranges
ranges = [0.3,0.2];
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



% constraint theta is processed: 2017-01-19 08:56:50
% constraint psy is processed: 2017-01-19 08:56:53
% constraint zeta is processed: 2017-01-19 08:56:55
% Optimization terminated.
% --------------------------------------------------------------
% The parameter setting:
% degree: 2; lambda degree: 0; eps1: 0.01; eps2: 0.01
% --------------------------------------------------------------
% The coefficients of function phy is:
%    -0.0856    0.1010   -0.7600    0.2408    0.7600   -0.0793
% 
% --------------------------------------------------------------
% The function phy is:
% (13875274223*x)/137438953472 - (13370672451591*y)/17592186044416 + (6685336225795*x*y)/8796093022208 + (2117713315515*x^2)/8796093022208 - (1394370808931*y^2)/17592186044416 - 1506710086745/17592186044416
%  
% --------------------------------------------------------------
% The coefficients of lambda is:
%      0
% 
% --------------------------------------------------------------
% The function lambda is:
% 0
%  
% --------------------------------------------------------------
% The rou is:
%    -0.0320
% 
% --------------------------------------------------------------
% The computation time is:
%     0.0710
% 
% --------------------------------------------------------------
% constraint theta is processed: 2017-01-19 08:57:03
% constraint psy is processed: 2017-01-19 08:57:06
% constraint zeta is processed: 2017-01-19 08:57:08
% Optimization terminated.
% Verification succeeded.
% The parameter setting:phy degree: 2; eps1: 0.01; eps2: 0.01
% --------------------------------------------------------------
% The coefficients of function phy is:
%    1.0e+03 *
% 
%    -0.9844    0.9947    0.2045    0.1125    1.0343    0.5413
% 
% --------------------------------------------------------------
% The function phy is:
% (494944190703259*x^2)/4398046511104 + (568621861850137*x*y)/549755813888 + (8749598905630757*x)/8796093022208 + (4761463631148951*y^2)/8796093022208 + (7196073473344305*y)/35184372088832 - 4329284543843047/4398046511104
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0155
% 
% --------------------------------------------------------------
% Verify feasible solution succeed, norms ;
%    1.0e-10 *
% 
%     0.1359    0.0076    0.0158
% 
% --------------------------------------------------------------
% 
% ans = 
% 
%   LinearProgram4_v3 with properties:
% 
%                           indvars: [1x2 sym]
%                            degree: 2
%                               phy: [1x1 sym]
%                               eps: [0.0100 0.0100]
%                                 f: [2x1 sym]
%                           decvars: [1x114 sym]
%                             exprs: [1x5 Constraint]
%                     phyPolynomial: [1x1 lp4util.SymbolicPolynomial]
%                     pLambdaDegree: 0
%                 pLambdaPolynomial: [1x1 lp4util.SymbolicPolynomial]
%     pLambdaNormalizedSymbolicVars: []
%                     wSymbolicVars: [6x1 sym]
%                       wExpression: [1x1 sym]
%                            rouVar: [1x1 sym]
%                       pPartitions: [1024x1 lp4util.Partition]
%                 pLambdaPartitions: [1024x1 lp4util.Partition]
%                          linprogF: [1x114 double]
%                             theta: [1x2 sym]
%                               psy: [1x2 sym]
%                              zeta: [1x2 sym]
%                            rouInc: 0
% 
