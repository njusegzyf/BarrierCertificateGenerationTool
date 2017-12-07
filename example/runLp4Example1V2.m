function [lp, solveRes] = runLp4Example1V2() 

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
f = [-(11/2)*x2+x2^2;
    6*x1-x1^2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 2*x1-8;
theta2 = x2;
theta3 = x2-1;

g_theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = (x1-1)/4;
psy2 = (x2-1)/4;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1-1;
zeta2 = x2-2;

g_zeta = [zeta1, zeta2];

import lp4.LinearProgram4_v2
lp = LinearProgram4_v2(vars);

lp.f = f;
lp.eps = eps;

% Set the degree of phy
degree = 3;
pLambdaDegree = 1;
lp = lp.setDegreeAndInit(degree, pLambdaDegree);

% Note:
% degree = 4 and pLambdaDegree = 1 is OK.
% degree = 3 and pLambdaDegree = 1 is OK.

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


warning('on')

echo off;
end
