function [lp, solveRes, lpVer, solveResVer] = runLp4Example1FullV3() 

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

% Set the degree of phy and lambda
degree = 3;
pLambdaDegree = 1;

% Note:
% degree = 4 and pLambdaDegree = 1 is OK.
% degree = 3 and pLambdaDegree = 1 is OK.

import lp4util.Partition
phyRange = Partition(-0.5, 0.5);
pLambdaRange = Partition(-0.5, 0.5);
phyRangeInVerify = Partition(-10, 10);



import lp4.runAndVerifyWithLambdaV3
[lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithLambdaV3(...
vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree,...
phyRange, pLambdaRange, phyRangeInVerify);

warning('on')

echo off;
end
