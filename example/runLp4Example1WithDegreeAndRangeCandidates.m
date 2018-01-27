function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4Example1WithDegreeAndRangeCandidates()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example1Problem();



% set the degree of phy and lambda
degrees = [3];
pLambdaDegrees = [1, 2, 3];

% set the ranges
ranges = [1, 0.5, 0.3, 0.15, 0.1];
pRanges = repmat(10, 1, length(ranges)); % 这里 10 表示验证时 phy 的范围限制在 [-10, 10]
import lp4util.createRangeCandidates
[phyRanges, pLambdaRanges, phyRangesInVerify] = createRangeCandidates(ranges, ranges, pRanges);

% Note; For degrees = 3, pLambdaDegrees = 1, ranges = 1, Ok.



% run and verify
import lp4.runAndVerifyWithDegreeAndRangeCandidates
[lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithDegreeAndRangeCandidates(...
    vars, f, eps, g_theta, g_psy, g_zeta,...
    degrees, pLambdaDegrees,...
    phyRanges, pLambdaRanges, phyRangesInVerify);

warning('on')

echo off;
end



% Optimal solution found.
% 
% --------------------------------------------------------------
% The parameter setting:
% degree: 3; lambda degree: 1; eps1: 0.001; eps2: 0.001
% --------------------------------------------------------------
% The coefficients of function phy is:
%          0   -0.5454         0   -0.5454         0         0    0.3606   -0.1653         0   -0.4358
% 
% --------------------------------------------------------------
% The function phy is:
% (203022641163783*x1^3)/562949953421312 - (5956382437273943*x1^2*x2)/36028797018963968 - (2456075546844339*x1^2)/4503599627370496 - (2456075546844339*x1)/4503599627370496 - (1962834853899003*x2^3)/4503599627370496
%  
% --------------------------------------------------------------
% The coefficients of lambda is:
%      0     0     0
% 
% --------------------------------------------------------------
% The function lambda is:
% 0
%  
% --------------------------------------------------------------
% The rou is:
%    -0.4546
% 
% --------------------------------------------------------------
% The computation time is:
%     0.0235
% 
% --------------------------------------------------------------
% constraint theta is processed: 2018-01-08 21:14:27
% constraint psy is processed: 2018-01-08 21:14:36
% constraint zeta is processed: 2018-01-08 21:14:41
% 
% Optimal solution found.
% 
% Verification succeeded.
% The parameter setting:phy degree: 3; eps1: 0.001; eps2: 0.001
% --------------------------------------------------------------
% The coefficients of function phy is:
%    -0.0115   -0.0002    0.0001    0.0013         0    0.0011   -0.0001         0         0   -0.0001
% 
% --------------------------------------------------------------
% The function phy is:
% - (5221541966922293*x1^3)/36893488147419103232 + (5874234712787581*x1^2)/4611686018427387904 - x1/4500 - (2337485886072857*x2^3)/18446744073709551616 + (5006193254385075*x2^2)/4611686018427387904 + (8833279599453753*x2)/73786976294838206464 - 6626057296624775/576460752303423488
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0130
% 
% --------------------------------------------------------------
% Verify feasible solution succeed, norms ;
%    1.0e-17 *
% 
%     0.3917    0.7744    0.1278
% 
% --------------------------------------------------------------
% 
% ans = 
% 
%   LinearProgram4_v3 - 属性:
% 
%                           indvars: [1×2 sym]
%                            degree: 3
%                               phy: [1×1 sym]
%                               eps: [1.0000e-03 1.0000e-03]
%                                 f: [2×1 sym]
%                           decvars: [1×233 sym]
%                             exprs: [1×5 Constraint]
%                     phyPolynomial: [1×1 lp4util.SymbolicPolynomial]
%                     pLambdaDegree: 1
%                 pLambdaPolynomial: [1×1 lp4util.SymbolicPolynomial]
%     pLambdaNormalizedSymbolicVars: []
%                     wSymbolicVars: [10×3 sym]
%                       wExpression: [1×1 sym]
%                            rouVar: [1×1 sym]
%                       pPartitions: [1024×1 lp4util.Partition]
%                 pLambdaPartitions: [1024×1 lp4util.Partition]
%                          linprogF: [1×233 double]
%                             theta: [1×3 sym]
%                               psy: [1×2 sym]
%                              zeta: [1×2 sym]
%                            rouInc: 0
