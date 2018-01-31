function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4Example5WithDegreeAndRangeCandidates()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example5Problem();
                              

% set the degree of phy and lambda
degrees = [2];
pLambdaDegrees = [0,1,2];

% set the ranges
ranges = [0.9,0.7,0.5,0.3,0.1];
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


% constraint theta is processed: 2017-12-20 15:08:21
% constraint psy is processed: 2017-12-20 15:08:27
% constraint zeta is processed: 2017-12-20 15:08:31
% Optimization terminated.
% --------------------------------------------------------------
% The parameter setting:
% degree: 3; lambda degree: 0; eps1: 0.001; eps2: 0.001
% --------------------------------------------------------------
% The coefficients of function phy is:
%     0.0982    0.7761   -0.5046    0.1927   -0.7761   -0.3652    0.0371   -0.7761   -0.7761   -0.0000
% 
% --------------------------------------------------------------
% The function phy is:
% (3413121114257*x1)/4398046511104 - (138707165353*x2)/274877906944 - (1706560557045*x1*x2)/2199023255552 - (3413121114675*x1*x2^2)/4398046511104 - (1706560557387*x1^2*x2)/2199023255552 + (211844392069*x1^2)/1099511627776 + (81531792741*x1^3)/2199023255552 - (1605981630487*x2^2)/4398046511104 - (619*x2^3)/4398046511104 + 215999294473/2199023255552
%  
% --------------------------------------------------------------
% The coefficients of lambda is:
%    2.9331e-11
% 
% --------------------------------------------------------------
% The function lambda is:
% 129/4398046511104
%  
% --------------------------------------------------------------
% The rou is:
%    -0.0192
% 
% --------------------------------------------------------------
% The computation time is:
%     0.2963
% 
% --------------------------------------------------------------
% constraint theta is processed: 2017-12-20 15:08:51
% constraint psy is processed: 2017-12-20 15:08:56
% constraint zeta is processed: 2017-12-20 15:09:00
% Optimization terminated.
% Verification succeeded.
% The parameter setting:phy degree: 3; eps1: 0.001; eps2: 0.001
% --------------------------------------------------------------
% The coefficients of function phy is:
%    1.0e+04 *
% 
%     0.1277   -0.0050    0.0847    0.0374   -1.1544    0.0559   -0.0140    0.4084    0.3967   -0.0130
% 
% --------------------------------------------------------------
% The function phy is:
% - (2458318146316985*x1^3)/17592186044416 + (4490175695470437*x1^2*x2)/1099511627776 + (3286747227913345*x1^2)/8796093022208 + (2180985540963769*x1*x2^2)/549755813888 - (3173061827940337*x1*x2)/274877906944 - (440805590092881*x1)/8796093022208 - (2282638084648505*x2^3)/17592186044416 + (4914647779357659*x2^2)/8796093022208 + (930760800380479*x2)/1099511627776 + 5617403996566297/4398046511104
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0206
% 
% --------------------------------------------------------------
% Verify feasible solution succeed, norms ;
%    1.0e-10 *
% 
%     0.1836    0.0679    0.1491
% 
% --------------------------------------------------------------
% 
% ans = 
% 
%   LinearProgram4_v3 with properties:
% 
%                           indvars: [1x2 sym]
%                            degree: 3
%                               phy: [1x1 sym]
%                               eps: [1.0000e-03 1.0000e-03]
%                                 f: [2x1 sym]
%                           decvars: [1x162 sym]
%                             exprs: [1x5 Constraint]
%                     phyPolynomial: [1x1 lp4util.SymbolicPolynomial]
%                     pLambdaDegree: 0
%                 pLambdaPolynomial: [1x1 lp4util.SymbolicPolynomial]
%     pLambdaNormalizedSymbolicVars: []
%                     wSymbolicVars: [10x1 sym]
%                       wExpression: [1x1 sym]
%                            rouVar: [1x1 sym]
%                       pPartitions: [1024x1 lp4util.Partition]
%                 pLambdaPartitions: [1024x1 lp4util.Partition]
%                          linprogF: [1x162 double]
%                             theta: [1x2 sym]
%                               psy: [1x2 sym]
%                              zeta: [1x2 sym]
%                            rouInc: 0
% 
