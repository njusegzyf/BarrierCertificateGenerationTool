function [lp, solveRes] = runLp4Example6()

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% independent variables
syms v a d vf af;
vars = [v a d vf af];

% Constructing the vector field dx/dt = f
f = [a;
     -3*a-3*(v-vf)+(d-(v+10));
     vf-v;
     af;
     0];

eps = [0.001,0.001];

% Constructing the theta constraint
theta1 = d-4;
theta2 = d-5;
theta3 = v-vf+1;
theta4 = v-vf;
theta5 = a+1;
theta6 = a;
g_theta = [theta1,theta2,theta3,theta4,theta5,theta6];

% Constructing the psy constraint
psy1 = d;
psy2 = v;
psy3 = vf;
g_psy = [psy1, psy2, psy3];

% Constructing the zeta constraint
zeta1 = (a+2)/7;
zeta2 = (af+2)/7;
zeta3 = d+1;
zeta4 = d;
g_zeta = [zeta1,zeta2,zeta3,zeta4];

import lp4.LinearProgram4_v2
lp = LinearProgram4_v2(vars);

lp.f = f;
lp.eps = eps;

% Set the degree of phy
degree = 3;
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
