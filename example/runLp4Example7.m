function [lp, solveRes] = runLp4Example7()

% Example 1.5 

clear; 
echo on;

% independent variables
syms x y;
vars = [x y];

% Constructing the vector field dx/dt = f
f = [1+x^2*y-2.5*x;
     1.5*x-x^2*y];

% Set the degree of phy
degree = 2;

M = 0;
r = 0;
eps = [0.1,0.1];

% Constructing the theta constraint
theta1 = 10*x-9;
theta2 = 10*y;

g_theta = [theta1,theta2];

% Constructing the psy constraint
psy1 = x;
psy2 = y;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 5*x-1;
zeta2 = 5*y-2;

g_zeta = [zeta1,zeta2];

import lp4.LinearProgram4
lp = LinearProgram4(vars);

lp.f = f;
lp.r = r;
lp.eps = eps;

lp = lp.setDegreeAndInit(degree);

lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

import lp4util.Partition
lp.pPartitions = repmat(Partition(-1200, 1200), 15, 1);
lp.pLambdaPartitions = repmat(Partition(-10, 10), 15, 1);
lp = lp.setWConstraint();

lp = lp.setDevVarsConstraint();

lp = lp.setLinprogF();

% solve the lp problem
[lp, solveRes] = lp.solve();

% constraint theta is generated: 2015-12-24 11:20:29
% constraint psy is generated: 2015-12-24 11:20:33
% constraint z  eta is generated: 2015-12-24 11:20:34
% constraint theta is processed: 2015-12-24 11:20:36
% constraint psy is processed: 2015-12-24 11:20:48
% constraint zeta is processed: 2015-12-24 11:20:51
% Optimization terminated.
% Elapsed time is 0.010418 seconds.
% --------------------------------------------------------------
% The parameter setting:
% degree: 2; r: 0; eps1: 0; eps1: 0.1; eps2: 0.1
% --------------------------------------------------------------
% The coefficients of function phy is:
%    1.0e+03 *
% 
%    -0.9856
%     0.9956
%     0.2042
%     0.1130
%     1.0366
%     0.5424
% 
% --------------------------------------------------------------
% The function phy is:
% (7955026750979167*x^2)/70368744177664 + (1139706653757397*x*y)/1099511627776 + (8757534132573727*x)/8796093022208 + (4770616731947699*y^2)/8796093022208 + (1796018599446961*y)/8796093022208 - 4334693643668357/4398046511104
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0105
% 
% --------------------------------------------------------------


echo off;
end