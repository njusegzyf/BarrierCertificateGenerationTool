function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4BenkEx17WithDegreeAndRangeCandidates()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx17Problem();

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [-1 + x^2 + y^2;
     5 * (-1 + x*y)];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = -x;
theta2 = -y;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 2) / 4;
psy2 = (y + 2) / 4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x - 1;
zeta2 = y - 1;
g_zeta = [zeta1, zeta2];



% Set the degree of phy and lambda
degrees = [1, 2, 3, 4];
pLambdaDegrees = [1, 2, 3];

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
