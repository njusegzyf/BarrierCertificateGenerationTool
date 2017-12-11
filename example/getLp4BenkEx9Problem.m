function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx9Problem()

% dai

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [2 * x - x * y;
     2 * x^2 - y];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = (x + 1) / 2;
theta2 = (y + 3) / 2;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 3) / 6;
psy2 = (y + 3) / 6;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x + 0.3) / 0.6;
zeta2 = (y - 0.7) / 0.6;
g_zeta = [zeta1, zeta2];

end
