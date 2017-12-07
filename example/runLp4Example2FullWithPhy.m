function [lp, solveRes, lpVer, solveResVer] = runLp4Example2FullWithPhy() 

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

% Set the degree of phy and lambda
degree = 3;
pLambdaDegree = 1;



import lp4.LinearProgram4_v2
lp = LinearProgram4_v2(vars);

lp.f = f;
lp.eps = eps;

lp = lp.setDegreeAndInit(degree, pLambdaDegree);

lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

import lp4util.Partition
lp.pPartitions = repmat(Partition(-1, 1), 1024, 1);
lp.pLambdaPartitions = repmat(Partition(-1, 1), 1024, 1);
lp = lp.setWConstraint();

lp = lp.setDevVarsConstraint();

lp = lp.setLinprogF();

% solve the lp problem
[lp, solveRes] = lp.solve();
% [lp, solveRes] = lp.solve1And3();

if ~(solveRes.hasSolution())
    disp('No solution to verify.');
    lpVer = 0;
    solveResVer = 0;
    return;
end
    
% verify

import lp4.LinearProgram4Verification3
lpVer = LinearProgram4Verification3(vars);

lpVer.f = f;
lpVer.eps = eps;

% Set the degree of phy
lpVer = lpVer.setDegreeAndInit(pLambdaDegree);

lpVer.phy = solveRes.getPhyExpression();

lpVer = lpVer.setThetaConstraint(g_theta);
lpVer = lpVer.setPsyConstraint(g_psy);
lpVer = lpVer.setZetaConstraint(g_zeta);
lpVer = lpVer.generateEqsForConstraint1To3();

lpVer = lpVer.setDevVarsConstraint();

[lpVer, solveResVer, resNorms] = lpVer.solve();

if ~(solveResVer.hasSolution())
    disp('Verify failed.');
else
    disp('Verify succeed, norms ;');
    disp(resNorms);
end

warning('on')

echo off;
end
