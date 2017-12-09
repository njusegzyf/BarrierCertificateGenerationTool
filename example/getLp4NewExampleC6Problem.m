function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC4Problem()

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [x2;
    -x1+(1/3)*x1^3-x2];

eps = [0.00001,0.00001];

% Constructing the theta constraint
theta1 = x1 - 1;
theta2 = x2 + 0.5;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x1 + 4) / 8;
psy2 = (x2 + 4) / 8;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 1.25 * x1 + 1.75;
zeta2 = 1.25 * x2 + 1.75;
g_zeta = [zeta1, zeta2];

end
