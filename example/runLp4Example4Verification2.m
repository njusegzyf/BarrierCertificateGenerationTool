function [lp, solveRes, verified]  = runLp4Example4Verification2()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [x1^2-2*x1+x2;
     x1+x2^2-2*x2;];

eps = [0.1,0.1];

% Constructing the theta constraint
theta1 = 100*(x1^2+x2^2);

g_theta = [theta1];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x1^2+x2^2-0.25)/7.75;

g_zeta = [zeta1];

import lp4.LinearProgram4Verification2
lp = LinearProgram4Verification2(vars);

lp.f = f;
lp.eps = eps;

% Note; degree >= 2

% Set the degree of phy
lp = lp.setDegreeAndInit(4);

lp.lambda = 0;
% - x2^2/2 - 1/2;

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
