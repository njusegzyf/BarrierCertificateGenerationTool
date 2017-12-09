function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC5Problem()

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [-x1 + 2 * x1^2 * x2;
     -x2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = x1 + 0.5;
theta2 = x1 - 0.5;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x1 + 2) / 4;
psy2 = (x2 + 2) / 4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1 - 1;
zeta2 = x2 - 1;
g_zeta = [zeta1, zeta2];

end
