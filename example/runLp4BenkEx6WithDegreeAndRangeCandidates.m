function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4BenkEx6WithDegreeAndRangeCandidates()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx6Problem();



% Set the degree of phy and lambda
degrees = [2, 3, 4];
pLambdaDegrees = [0, 1, 2];

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
% degree: 2; lambda degree: 0; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function phy is:
%     0.1868         0    0.4729   -0.0000         0   -0.0390
% 
% --------------------------------------------------------------
% The function phy is:
% (8519145237142461*y)/18014398509481984 - (894765698930641*x^2)/36893488147419103232 - (2812262120919233*y^2)/72057594037927936 + 6728396027305761/36028797018963968
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
%    -0.0542
% 
% --------------------------------------------------------------
% The computation time is:
%     0.0251
% 
% --------------------------------------------------------------
% constraint theta is processed: 2017-12-11 12:22:51
% constraint psy is processed: 2017-12-11 12:23:11
% constraint zeta is processed: 2017-12-11 12:23:15
% 
% No feasible solution found.
% 
% Linprog stopped because no point satisfies the constraints.
% 
% --------------------------------------------------------------
% The parameter setting:
% degree: 2; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The problem with degree 2 maybe have no solution.
% --------------------------------------------------------------
% Verify feasible solution failed.
% constraint theta is processed: 2017-12-11 12:23:26
% constraint psy is processed: 2017-12-11 12:23:45
% constraint zeta is processed: 2017-12-11 12:23:49
% 
% Optimal solution found.
% 
% --------------------------------------------------------------
% The parameter setting:
% lambdaDegree: 0; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function lambda is:
%    -0.7769
% 
% --------------------------------------------------------------
% The function lambda is:
% -6997937753528645/9007199254740992
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0177
% 
% --------------------------------------------------------------
% Verify feasible solution succeed, norms ;
%    1.0e-15 *
% 
%     0.0621    0.0636    0.1146
