function [lp, solveRes, verified] = runLp4Example3Verification2()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [2*x1-x1*x2;
     2*x1^2-x2;];

eps = [0.001,0.001];

% Constructing the theta constraint
theta1 = x1^2+(x2+2)^2;

g_theta = [theta1];

% Constructing the psy constraint
psy1 = (x1+7)/14;
psy2 = (x2+7)/14;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1^2+(x2-5.2)^2;

g_zeta = [zeta1];

import lp4.LinearProgram4Verification2
lp = LinearProgram4Verification2(vars);

lp.f = f;
lp.eps = eps;

% Set the degree of phy
lp = lp.setDegreeAndInit(3);

lp.lambda = - (26617*x1)/81520 - (26617*x2)/81520 - 10211/20380;

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
