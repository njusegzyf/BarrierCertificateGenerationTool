function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftTSC_C8Problem()

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [-x1 + 2 * x1^2 * x2;
     -x2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = (x1 + 1/4) * 2;
theta2 = (x2 - 3/4) * 4/3;
g_theta = [theta1, theta2];

% Constructing the psy constraint
g_psy = (vars + 2)/4;

% Constructing the zeta constraint
g_zeta = vars - 1;

end
