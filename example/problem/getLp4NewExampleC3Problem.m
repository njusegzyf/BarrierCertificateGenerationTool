function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC3Problem()

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [x1 + x2;
     x1 * x2 - 0.5 * x2^2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 5 * x1;
theta2 = 5 * x2 + 1;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x1 + 5) * (4/21);
psy2 = (x2 + 5) * 0.1;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
sqrtV = sqrt(0.34);
zeta1 = (x1 + 1.5 + sqrtV) / (2 * sqrtV);
zeta2 = (x2 + 1.5 + sqrtV) / (2 * sqrtV);
g_zeta = [zeta1, zeta2];

end
