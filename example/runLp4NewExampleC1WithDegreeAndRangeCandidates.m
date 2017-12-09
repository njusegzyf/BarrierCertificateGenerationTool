function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4NewExampleC1WithDegreeAndRangeCandidates()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC1Problem();

% Set the degree of phy and lambda
degrees = [1, 2, 3, 4];
pLambdaDegrees = [1, 2, 3];

rangeCandidates = [1, 0.5, 0.3, 0.15, 0.1];
import lp4util.Partition
phyRanges = arrayfun(@(x) Partition(-x, x), rangeCandidates);
pLambdaRanges = arrayfun(@(x) Partition(-x, x), rangeCandidates);
phyRangesInVerify = 0;

% Note: For degree = [], pLambdaDegree = [], range = [1..0.1], unable to verify feasible solutions.

import lp4.runAndVerifyWithDegreeAndRangeCandidates
[lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithDegreeAndRangeCandidates(...
    vars, f, eps, g_theta, g_psy, g_zeta,...
    degrees, pLambdaDegrees,...
    phyRanges, pLambdaRanges, phyRangesInVerify);

warning('on')

echo off;
end
