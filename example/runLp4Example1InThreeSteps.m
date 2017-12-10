function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runLp4Example1InThreeSteps() 

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
theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = (x1-1)/4;
psy2 = (x2-1)/4;
psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1-1;
zeta2 = x2-2;
zeta = [zeta1, zeta2];

% Set the degree of phy and lambda
degree = 3;
pLambdaDegree = 1;

% Note:
% degree = 4 and pLambdaDegree = 1 is OK.
% degree = 3 and pLambdaDegree = 1 is OK.

import lp4util.Partition
phyRange = Partition(-1, 1);
pLambdaRange = Partition(-1, 1);
phyRangeInVerify = 0; % Partition(-10, 10);

phyDelta = 0.1;
pLambdaDelta = 0.1;



import lp4.LinearProgram4_v3
lp = LinearProgram4_v3.createLp(vars, f, eps, theta, psy, zeta, degree, pLambdaDegree, phyRange, pLambdaRange);

% solve the lp problem
[lp, solveRes] = lp.solve();

if ~solveRes.hasSolution()
    return;
end

import lp4.runAndverifyFromSolveResult
[lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndverifyFromSolveResult(...
    lp, solveRes, phyDelta, pLambdaDelta, phyRangeInVerify);

warning('on')

echo off;
end
