function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC11Problem()

% independent variables
syms x1 x2 x3 x4;
vars = [x1, x2, x3, x4];

% Constructing the vector field dx/dt = f
f = [-x1 + x2^3 - 3 * x3 * x4;
     -x1 - x2^3;
	 x1 * x4 - x3;
	 x1 * x3 - x4^3];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 2 * x1;
theta2 = 2 * x2;
theta3 = 2 * x3;
theta4 = 2 * x4;
g_theta = [theta1, theta2, theta3, theta4];

% Constructing the psy constraint
psy1 = 0.5 * x1 + 0.5;
psy2 = 0.5 * x2 + 0.5;
psy3 = 0.5 * x3 + 0.5;
psy4 = 0.5 * x4 + 0.5;
g_psy = [psy1, psy2, psy3, psy4];

% Constructing the zeta constraint
zeta1 = 2 * x1 + 2;
zeta2 = 2 * x2 + 2;
zeta3 = 2 * x3 + 2;
zeta4 = 2 * x4 + 2;
g_zeta = [zeta1, zeta2, zeta3, zeta4];

end
