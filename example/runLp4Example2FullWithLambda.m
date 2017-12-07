function [lp, solveRes, lpVer, solveResVer] = runLp4Example2FullWithLambda() 

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
    -x1+(1/3)*x1^3-x2];

eps = [0.001,0.001];

% phy = 151/99+152/99*x1+62/33*x2+106/99*x1*x2+4/9*x1^2
% Constructing the theta constraint
theta1 = 4*(x1-1.5)^2+4*x2^2;
theta2 = x2+1/2;
theta3 = x1-1;

g_theta = [theta1,theta2,theta3];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 6.25*(x1+1)^2+6.25*(x2+1)^2;
zeta2 = x1+7/5;
zeta3 = x2+7/5;
g_zeta = [zeta1,zeta2,zeta3];

% Set the degree of phy and lambda
degree = 3;
pLambdaDegree = 7;

% Note: degree = 4 and pLambdaDegree = 1 is OK.

import lp4util.Partition
phyRange = Partition(-0.3, 0.3);
pLambdaRange = Partition(-0.3, 0.3);
phyRangeInVerify = 0; % Partition(-10, 10);



import lp4.runAndVerifyWithLambda
[lp, solveRes, lpVer, solveResVer, resNorms] = runAndVerifyWithLambda(...
vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree,...
phyRange, pLambdaRange, phyRangeInVerify);

warning('on')

echo off;
end
