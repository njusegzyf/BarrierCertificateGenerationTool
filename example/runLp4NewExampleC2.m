function [lp, solveRes, lpVer, solveResVer] = runLp4NewExampleC2()

% C2[19]

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [x1 - x2;
     x1 + x2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = x1 - 2.25;
theta2 = x2 + 0.5;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = 0.2 * x1;
psy2 = 0.2 * x2;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
sqrtV = sqrt(4.1);
zeta1 = (x1 - 1 + sqrtV) / (2 * sqrtV);
zeta2 = (x2 - 2 + sqrtV) / (2 * sqrtV);
g_zeta = [zeta1, zeta2];

% Set the degree of phy and lambda
degree = 3;
pLambdaDegree = 5;

import lp4util.Partition
phyRange = Partition(-1, 1);
pLambdaRange = Partition(-1, 1);
phyRangeInVerify = 0; % Partition(-10, 10);



import lp4.runAndVerifyWithLambdaV3
[lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithLambdaV3(...
vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree,...
phyRange, pLambdaRange, phyRangeInVerify);

warning('on')

echo off;
end
