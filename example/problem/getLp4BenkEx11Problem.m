function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx11Problem()

% darbowx

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [3 * x^2 - 12;
     x*y - y^2 + 3];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = 0.1 * x + 0.5;
theta2 = 0.1 * y + 0.5;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 3) / 6;
psy2 = (y + 3) / 6;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x + 5;
zeta2 = 0.1 * y + 0.5;
g_zeta = [zeta1, zeta2];

end
