function [lp, solveRes] = runLp4Example3()

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

import lp4.LinearProgram4_v2
lp = LinearProgram4_v2(vars);

lp.f = f;
lp.eps = eps;

% Set the degree of phy
degree = 2;
pLambdaDegree = 2;
lp = lp.setDegreeAndInit(degree, pLambdaDegree);



lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

import lp4util.Partition
lp.pPartitions = repmat(Partition(-1, 1), 1000, 1);
lp.pLambdaPartitions = repmat(Partition(-1, 1), 1000, 1);
lp = lp.setWConstraint();

lp = lp.setDevVarsConstraint();

lp = lp.setLinprogF();

% see http://blog.sina.com.cn/s/blog_68b0c65f0100mq5m.html
% options = optimset(¡®MaxIter¡¯, 2000);

% solve the lp problem
[lp, solveRes] = lp.solve();
% [lp, solveRes] = lp.solve(options);
% [lp, solveRes] = lp.solve1And3();

warning('on')

echo off;
end
