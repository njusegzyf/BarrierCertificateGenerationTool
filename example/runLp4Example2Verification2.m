function [lp, solveRes, verified] = runLp4Example2Verification2()

% Example 2

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

eps = [0.00001,0.00001];

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

import lp4.LinearProgram4Verification2
lp = LinearProgram4Verification2(vars);

lp.f = f;
lp.eps = eps;

% Set the degree of phy
lp = lp.setDegreeAndInit(4);

lp.lambda = - (2452868740955843*x1)/4503599627370496 - (919825777858441*x2)/1125899906842624 - (4905737481911685*x1*x2)/9007199254740992 - (7358606222867529*x1*x2^2)/9007199254740992 - (254984226781435*x1^2*x2)/281474976710656 - (462838824154013*x1^2)/562949953421312 - (5792728612319129*x1^3)/9007199254740992 - (4905737481911687*x2^2)/9007199254740992 - (1226434370477921*x2^3)/2251799813685248 - 2452868740955843/4503599627370496;
 
lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

lp = lp.setDevVarsConstraint();

% solve the lp problem
[lp, solveRes] = lp.solve();

if solveRes.hasSolution()
    verified = solveRes.verify();
else
    verified = false;
end

warning('on')

echo off;
end
