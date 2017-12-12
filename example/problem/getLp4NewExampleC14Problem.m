function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC14Problem()

% independent variables
syms x1 x2 x3;
vars = [x1, x2, x3];

% Constructing the vector field dx/dt = f
f = [x1^2 + x1 * x2 - x1 * x4;
     2 * x1 * x2 + x2^2;
	 x2 * x3 - 2 * x3^2];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = 2 * x1 - 2.5;
theta2 = 2 * x2 - 0.5;
theta3 = 2 * x3 - 2.5;
g_theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = 0.25 * x1 + 0.5;
psy2 = 0.25 * x2 + 0.5;
psy3 = 0.25 * x3 + 0.5;
g_psy = [psy1, psy2, psy3];

% Constructing the zeta constraint
zeta1 = 2 * x1 - 0.5;
zeta2 = 2 * x2 - 2.5;
zeta3 = 2 * x3 - 2.5;
g_zeta = [zeta1, zeta2, zeta3];

end
