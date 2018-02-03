function [vars, stateNum, fs, eps, thetaStateIndex, g_theta, g_psys, g_zetas, g_guards] = getLp4HybridEx3Problem() 

% Note: This problem is solved with degree = 1, degree = 1; pLambdaDegree = 1; pReDegree = 1;

% independent variables
syms x1 x2;
vars = [x1 x2];

stateNum = 2;

% Constructing the vector field dx/dt = f
f1 = [-0.1-0.2*x1; 1+0.2*x1-1.7*x2+0.8*x2^2];
f2 = [-0.3-0.2*x1+0.2*x2; 0.7-0.2*x1-0.2*x2];
fs = [f1, f2];
% Note: use `fs(:, 1)` to get f1

import lp4.Lp4Config
% eps = [Lp4Config.DEFAULT_EPS, Lp4Config.DEFAULT_EPS];
eps = [0.01, 0.01];


% Constructing the theta constraint
theta1 = (x1-5.25)/0.5;
theta2 = (x2)/0.5;
theta3 = ((x1-5.5)^2+(x2-0.25)^2)/0.0625;
g_theta = [theta1, theta2, theta3];
thetaStateIndex = 1;



% Constructing the psy constraint
psy11 = (x1-4)/2;
psy12 = x2;
g_psy1 = [psy11, psy12];

psy21 = (x1-4)/2;
psy22 = x2-1;
g_psy2 = [psy21, psy22];

g_psys = [g_psy1; g_psy2];



% Constructing the Guard constraint
import lp4.Guard

guard121=(x1-4)/2;
guard122=100*(x2-0.99);
guard12 = Guard(1, 2, [guard121, guard122], [x2], [1]);

g_guards = [guard12]; % [guard12, guard21, guard23, guard32];



% Constructing the zeta constraint
import lp4.UnsafeConstraint
zeta11 = ((x1-4.25)^2+(x2-0.25)^2)/0.0625;
zeta12 = (x1-4)/0.5;
zeta13 = x2/0.5;
g_zetas = [UnsafeConstraint(1, [zeta11, zeta12, zeta13])];


end
