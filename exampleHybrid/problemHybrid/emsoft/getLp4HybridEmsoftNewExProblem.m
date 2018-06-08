function [vars, stateNum, fs, eps, thetaStateIndex, g_theta, g_psys, g_zetas, g_guards] = getLp4HybridEmsoftNewExProblem()

% independent variables
syms x;
vars = [x];

stateNum = 2;

% Constructing the vector field dx/dt = f
fs = [-x, 30-x];
% Note: use `fs(:, 1)` to get f1

import lp4.Lp4Config
% eps = [Lp4Config.DEFAULT_EPS, Lp4Config.DEFAULT_EPS];
eps = [0.000001, 0.000001];
% eps = [0, 0];


% Constructing the theta constraint
theta1 = (x-19)/2;
g_theta = [theta1];
thetaStateIndex = 1;



% Constructing the psy constraint
psy11 = (x-18)/4;
g_psy1 = [psy11];

psy21 = (x-18)/4;
g_psy2 = [psy21];

g_psys = [g_psy1; g_psy2];



% Constructing the Guard constraint
import lp4.Guard

guard121=(x-18);
guard12 = Guard(1, 2, [guard121], [], []);

guard211=(x-21);
guard21 = Guard(2, 1, [guard211], [], []);

g_guards = [guard12, guard21];



% Constructing the zeta constraint
import lp4.UnsafeConstraint
zeta11 = (x-21.5)*2;
g_zetas = [UnsafeConstraint(2, [zeta11])];

end
