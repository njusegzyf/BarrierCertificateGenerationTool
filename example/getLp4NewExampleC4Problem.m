function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC4Problem()

% parillo from [23]

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [-x1 + x1 * x2;
     -x2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 4 * x1 - 4;
theta2 = 4 * x2 - 2;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x1 + 5) / 10;
psy2 = (x2 + 5) / 10;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x1 - 0.775) * 5;
zeta2 = (x2 - 0.05) * 5;
g_zeta = [zeta1, zeta2];

end
