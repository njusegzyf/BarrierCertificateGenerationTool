function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx10Problem()

% dai

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [y;
     2*x - x^3 - y - x^2 * y];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = (x + 1.4) / 0.8;
theta2 = (y - 1.6) / 0.8;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 3) / 6;
psy2 = (y + 3) / 6;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x - 0.8) / 0.4;
zeta2 = (y + 0.2) / 0.4;
g_zeta = [zeta1, zeta2];

end
