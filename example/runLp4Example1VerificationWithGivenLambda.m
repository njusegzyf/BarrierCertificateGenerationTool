function [lp, solveRes, resNorms]  = runLp4Example1VerificationWithGivenLambda()

% Example 1 

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

g_theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = (x1-1)/4;
psy2 = (x2-1)/4;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1-1;
zeta2 = x2-2;

g_zeta = [zeta1, zeta2];



phyDegree = 3;
lambda = - (1125925798617485*x1)/2251799813685248 - (1125925798617485*x2)/2251799813685248 - 1125925798617485/2251799813685248;



import lp4.LinearProgram4Verification2
lp = LinearProgram4Verification2(vars);

lp.f = f;
lp.eps = eps;

% Set the degree of phy
lp = lp.setDegreeAndInit(phyDegree);

lp.lambda = lambda;
% 0; % OK

lp = lp.setThetaConstraint(g_theta);
lp = lp.setPsyConstraint(g_psy);
lp = lp.setZetaConstraint(g_zeta);
lp = lp.generateEqsForConstraint1To3();

lp = lp.setDevVarsConstraint();

% solve the lp problem
[lp, solveRes, resNorms] = lp.solve();

if ~(solveRes.hasSolution())
    disp('Verify failed.');
else
    disp('Verify succeed, norms ;');
    disp(resNorms);
end

warning('on')

echo off;
end
