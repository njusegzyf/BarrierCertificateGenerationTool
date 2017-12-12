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
