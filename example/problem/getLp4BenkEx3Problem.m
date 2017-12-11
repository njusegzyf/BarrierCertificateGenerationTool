function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx3Problem()

% arrowswitch

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [x - y^3;
     x^3 + y];

eps = [0.000001, 0.000001];

% Constructing the theta constraint
theta1 = x - 1;
theta2 = y - 1;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 2) / 4;
psy2 = (y + 2) / 4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x + 0.5;
zeta2 = y + 0.5;
g_zeta = [zeta1, zeta2];

end
