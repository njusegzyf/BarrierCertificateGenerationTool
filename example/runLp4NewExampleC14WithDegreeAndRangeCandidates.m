function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4NewExampleC14WithDegreeAndRangeCandidates()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC14Problem();



% set the degree of phy and lambda
degrees = [2, 3, 4];
pLambdaDegrees = [0, 1];

% set the ranges
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
