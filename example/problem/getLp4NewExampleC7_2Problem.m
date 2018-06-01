function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4NewExampleC7_2Problem()

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [1+x1^2*x2-2.5*x1, 1.5*x1-x1^2*x2];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 10*(x1-0.9);
theta2 = 10*x2;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x1-0.2)/0.29;
zeta2 = (x2-0.3)/0.3;
g_zeta = [zeta1, zeta2];

end
