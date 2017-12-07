function [lp, solveRes] = runLp4NewExample1() 

% Example 1 

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [x2;
    -x1 + 1/3 * x1 ^ 3 - x2];

eps = [0.01, 0.01];

% Constructing the theta constraint
theta1 = x1-1;
theta2 = x2 + 0.5;

g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = 0.25 * x1 + 0.5;
psy2 = 0.25 * x2 + 0.5;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 1.25 * x1 + 1.75;
zeta2 = 1.25 * x2 + 1.75;

g_zeta = [zeta1, zeta2];

import lp4.LinearProgram4
lp = LinearProgram4(vars);

lp.f = f;
lp.eps = eps;

% Set the degree of phy
degree = 4;
pLambdaDegree = 1;
lp = lp.setDegreeAndInit(degree, pLambdaDegree);

% Note:
% degree = 4 and pLambdaDegree = 1 is OK for condition 1 and 3.

lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

import lp4util.Partition
lp.pPartitions = repmat(Partition(-1, 1), 15, 1);
lp.pLambdaPartitions = repmat(Partition(-1, 1), 15, 1);
lp = lp.setWConstraint();

lp = lp.setDevVarsConstraint();

lp = lp.setLinprogF();

% solve the lp problem
[lp, solveRes] = lp.solve();

% [lp, solveRes] = lp.solve1And3();

warning('on')

echo off;
end
