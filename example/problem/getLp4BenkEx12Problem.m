function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx11Problem()

% darbowx

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [x + y;
     x * y - y^2 / 2];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = (x + 1.45) / 0.4;
theta2 = (y - 1.05) / 0.4;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 3) / 6;
psy2 = (y + 3) / 6;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x + 2.7) / 0.4;
zeta2 = (y - 0.6) / 0.4;
g_zeta = [zeta1, zeta2];

end
