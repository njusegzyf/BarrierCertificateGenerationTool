function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4NewExampleC1WithRangeCandidates()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [-x1^3 / 3 + x1 - x2 + 7/8;
     2/25 * (x1 - 0.8 * x2 + 0.7)];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 2 * x1 + 2;
theta2 = 2 * x2 - 2;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = x1/8 + 0.5;
psy2 = x2/8 + 0.5;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 2 * x1 + 5;
zeta2 = 2 * x2 + 4;
g_zeta = [zeta1, zeta2];

% Set the degree of phy and lambda
degree = 2;
pLambdaDegree = 2;

rangeCandidates = [1, 0.5, 0.3, 0.15, 0.1];
import lp4util.Partition
phyRanges = arrayfun(@(x) Partition(-x, x), rangeCandidates);
pLambdaRanges = arrayfun(@(x) Partition(-x, x), rangeCandidates);
phyRangesInVerify = 0;

% Note: For degree = 2, pLambdaDegree = 1, range from 1 to 0.1, unable to verify feasible solutions.


import lp4.runAndVerifyWithRangeCandidates
[lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithRangeCandidates(...
vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree,...
phyRanges, pLambdaRanges, phyRangesInVerify);

warning('on')

echo off;
end
