function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC16Problem()

% independent variables
syms x1 x2 x3 x4;
vars = [x1, x2, x3, x4];

% Constructing the vector field dx/dt = f
f = [-0.5 * x1^2 - 0.5 * x2^2 - 0.125 * x3^2 - 2 * x2 * x3 + 2 * x4^2 + 1;
     -x1 * x2 - 1;
	 -x1 * x3;
	 -x1 * x4];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = 2 * x1 - 1.5;
theta2 = 2 * x2 - 1.5;
theta3 = 2 * x3 - 1.5;
theta4 = 2 * x4 - 1.5;
g_theta = [theta1, theta2, theta3, theta4];

% Constructing the psy constraint
psy1 = 0.25 * x1 + 0.5;
psy2 = 0.25 * x2 + 0.5;
psy3 = 0.25 * x3 + 0.5;
psy4 = 0.25 * x4 + 0.5;
g_psy = [psy1, psy2, psy3, psy4];

% Constructing the zeta constraint
zeta1 = 2 * x1 - 2.5;
zeta2 = 2 * x2 - 2.5;
zeta3 = 2 * x3 - 2.5;
zeta4 = 2 * x4 - 2.5;
g_zeta = [zeta1, zeta2, zeta3, zeta4];

end
