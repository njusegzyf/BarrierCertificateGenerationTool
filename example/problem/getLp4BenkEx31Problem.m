function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx29Problem()

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [-2 * y;
     -2 * x - 3 * x^2];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = 2.5 * (x + 13/15);
theta2 = (y + 0.2) / 0.4;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 2) / 4;
psy2 = (y + 2) / 4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x - 1;
zeta2 = y - 1;
g_zeta = [zeta1, zeta2];

end
