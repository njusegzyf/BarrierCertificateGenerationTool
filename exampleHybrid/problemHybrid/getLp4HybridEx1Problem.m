function [vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx1Problem()

% independent variables
syms T t;
vars = [T, t];
% vars = [t, T]; % wrong

stateNum = 3;

% Constructing the vector field dx/dt = f
f1 = [-T; 1];  % cool
f2 = [2; 1];    %heat
f3 = [-T/2; 1];  %check
fs = [f1, f2, f3];
% Note: use `fs(:, 1)` to get f1

import lp4.Lp4Config
% eps = [Lp4Config.DEFAULT_EPS, Lp4Config.DEFAULT_EPS];
eps = [0.1, 0.1];


% Constructing the theta constraint
theta1 = t;
theta2 = -t;
theta3 = (T-5)/5;
theta = [theta1, theta2, theta3];
thetaStateIndex = 2;



% Constructing the psy constraint
psy11 = t/100;
psy12 = (T-5)/95;
g_psy1 = [psy11, psy12];

psy21 = t/3;
psy22 = T/10;
g_psy2 = [psy21, psy22];

psy31 = t;
psy32 = T/100;
g_psy3 = [psy31, psy32];

psys = [g_psy1; g_psy2; g_psy3];



% Constructing the Guard constraint
import lp4.Guard

guard121=T/6;
guard122=t/100;
% `[t], [0]` means reset `t` to `0`
guard12 = Guard(1, 2, [guard121, guard122], [t], [0]);

guard211=(T-9)/91;
guard212=t/100;
guard21 = Guard(2, 1, [guard211, guard212], [], []);

guard231=T/100;
guard232=(t-2)/98;
guard23 = Guard(2, 3, [guard231, guard232], [t], [0]);

guard321=T/100;
guard322=(t-0.5)/99.5;
guard32 = Guard(3, 2, [guard321, guard322], [t], [0]);

guards = [guard12, guard21, guard23, guard32];



% Constructing the zeta constraint
import lp4.UnsafeConstraint
zeta31 = T/4.5;
zeta32 = t/100;
zetas = [UnsafeConstraint(3, [zeta31, zeta32])];

end
