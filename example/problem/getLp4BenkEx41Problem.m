function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx41Problem()

% wigging

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [-x + 2 * y + x^2 * y + x^4 * y^5;
     -y - x^4 * y^6 + x^8 * y^9];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = 2 * (x - 1/12);
theta2 = 2 * (y - 1/12);
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 2) / 4;
psy2 = (y + 2) / 4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x + 2;
zeta2 = y + 2;
g_zeta = [zeta1, zeta2];

end
