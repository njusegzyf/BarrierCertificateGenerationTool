function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx20Problem()

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [y;
     2 * (-1 - y) * y];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = -2 * x;
theta2 = -2 * y;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 3) / 6;
psy2 = (y + 3) / 6;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x - 2;
zeta2 = y - 2;
g_zeta = [zeta1, zeta2];

end
