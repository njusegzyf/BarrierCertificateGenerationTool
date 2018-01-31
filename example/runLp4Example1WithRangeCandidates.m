function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4Example1WithRangeCandidates()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example1Problem();

% Set the degree of phy and lambda
degree = 3;
pLambdaDegree = 0;

rangeCandidates = [1, 0.5, 0.3, 0.15, 0.1];
import lp4util.Partition
phyRanges = arrayfun(@(x) Partition(-x, x), rangeCandidates);
pLambdaRanges = arrayfun(@(x) Partition(-x, x), rangeCandidates);
phyRangesInVerify = [];

% Note: degree = 3, pLambdaDegree = 0 is Ok.



import lp4.runAndVerifyWithRangeCandidatesV2
[lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithRangeCandidatesV2(...
    vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree,...
    phyRanges, pLambdaRanges, phyRangesInVerify);

warning('on')

echo off;
end
