function [lp, solveRes, lpVer, solveResVer] = runLp4NewExample2()

% simdim from [16]

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% independent variables
syms x1 x2 x3 x4 x5 x6;
vars = [x1, x2, x3, x4, x5, x6];

% Constructing the vector field dx/dt = f
f = [-3 * x1^3 + 4 * x2^3 - 6 * x3 * x4;
    -x1 - x2 + x5^3;
    x1 * x4 - x3 + x4 * x6;
    x1 * x3 + x3 * x6 - x4^3;
    -2 * x2^3 - x5 + x6;
    -3 * x3 * x4 - x5^3 - x6];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 10 * x1 - 30;
theta2 = 10 * x2 - 30;
theta3 = 10 * x3 - 30;
theta4 = 10 * x4 - 30;
theta5 = 10 * x5 - 30;
theta6 = 10 * x6 - 30;
g_theta = [theta1, theta2, theta3, theta4, theta5, theta6];

% Constructing the psy constraint
psy1 = 10 * x1 - 40;
psy2 = 10 * x2 - 41;
psy3 = 10 * x3 - 42;
psy4 = 10 * x4 - 43;
psy5 = 10 * x5 - 44;
psy6 = 10 * x6 - 45;
g_psy = [psy1, psy2, psy3, psy4, psy5, psy6];

% Constructing the zeta constraint
zeta1 = x1 / 10;
zeta2 = x2 / 10;
zeta3 = (x3 - 2) / 8;
zeta4 = x4 / 10;
zeta5 = x5 / 10;
zeta6 = x6 / 10;
g_zeta = [zeta1, zeta2, zeta3, zeta4, zeta5, zeta6];

% Set the degree of phy and lambda
degree = 3;
pLambdaDegree = 1;

import lp4util.Partition
phyRange = Partition(-0.3, 0.3);
pLambdaRange = Partition(-0.3, 0.3);
phyRangeInVerify = 0; % Partition(-100, 100);



import lp4.runAndVerifyWithLambda
[lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithLambda(...
vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree,...
phyRange, pLambdaRange, phyRangeInVerify);

warning('on')

echo off;
end
