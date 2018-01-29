function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4HybridEx1State3WithCandidates()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4HybridEx1State3Problem();



% set the degree of phy and lambda
degrees =  [4, 5]; % [1, 2, 3, 4];
pLambdaDegrees = [1]; % [1, 2, 3];

% set the ranges
ranges = [1, 0.3, 0.1];
import lp4util.createRangeCandidates
[phyRanges, pLambdaRanges, ~] = createRangeCandidates(ranges, ranges, []);


% run and verify
import lp4.runAndVerifyWithDegreeAndRangeCandidates
[lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithDegreeAndRangeCandidates(...
    vars, f, eps, g_theta, g_psy, g_zeta,...
    degrees, pLambdaDegrees,...
    phyRanges, pLambdaRanges, []);

warning('on')

echo off;
end
