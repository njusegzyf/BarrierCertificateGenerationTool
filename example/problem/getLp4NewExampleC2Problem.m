function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC1Problem()

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [x1 - x2;
     x1 + x2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = x1 - 2.25;
theta2 = x2 + 0.5;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = 0.2 * x1;
psy2 = 0.2 * x2;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
% sqrtV = sqrt(4.1);
% zeta1 = (x1 - 1 + sqrtV) / (2 * sqrtV);
% zeta2 = (x2 - 2 + sqrtV) / (2 * sqrtV);
% Note: Unsafe set is changed to avoid intersecting with init set.
zeta1 = x1 / 2;
zeta2 = (x2 - 1) / 2;
g_zeta = [zeta1, zeta2];

end
