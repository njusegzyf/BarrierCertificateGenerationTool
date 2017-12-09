function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4NewExampleC1WithRangeCandidates()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC1Problem();

% Set the degree of phy and lambda
degree = 1;
pLambdaDegree = 2;

rangeCandidates = [1, 0.5, 0.3, 0.15, 0.1];
import lp4util.Partition
phyRanges = arrayfun(@(x) Partition(-x, x), rangeCandidates);
pLambdaRanges = arrayfun(@(x) Partition(-x, x), rangeCandidates);
phyRangesInVerify = 0;

% Note: For degree = 1, pLambdaDegree = 1, range from 1 to 0.1, unable to verify feasible solutions.
% Note: For degree = 2, pLambdaDegree = 1, range from 1 to 0.1, unable to verify feasible solutions.
% Note: For degree = 2, pLambdaDegree = 2, range from 1 to 0.1, unable to verify feasible solutions.
% Note: For degree = 3, pLambdaDegree = 1, range from 1 to 0.1, unable to verify feasible solutions.

import lp4.runAndVerifyWithRangeCandidates2
[lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithRangeCandidates2(...
vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree,...
phyRanges, pLambdaRanges, phyRangesInVerify);

warning('on')

echo off;
end
