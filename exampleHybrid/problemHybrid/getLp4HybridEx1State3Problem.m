function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4HybridEx1State3Problem()

% independent variables
syms T t;
vars = [T, t];

% Constructing the vector field dx/dt = f
f = [-T/2; 1];  %check

import lp4.Lp4Config
eps = [Lp4Config.DEFAULT_EPS, Lp4Config.DEFAULT_EPS];

% Constructing the theta constraint
theta1 = t;
theta2 = -t;
theta3 = (T - 5) / 5; % T / 10;
g_theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = t;
psy2 = T/100;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta31 = T/4.5;
zeta32 = t/100;
g_zeta = [zeta31, zeta32];

end
