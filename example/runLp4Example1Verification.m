function [lp, solveRes]  = runLp4Example1Verification()

% Example 1 

clear; 
echo on;

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [-(11/2)*x2+x2^2;
    6*x1-x1^2];

% Set the degree of phy
degree = 4;

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

import lp4.LinearProgram4Verification
lp = LinearProgram4Verification(vars);

lp.f = f;
lp.eps = eps;

lp.degree = degree;

lp.phy = (5376570161281769*x1^2)/2361183241434822606848 - (1351185711459749*x2)/73786976294838206464 - (1337099369181137*x2^2)/73786976294838206464;
lp.lambda = - x1^2 - x1*x2 - x1 - x2^2 - x2 - 1;

lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

lp = lp.setDevVarsConstraint();

% solve the lp problem
[lp, solveRes] = lp.solve();

echo off;
end
