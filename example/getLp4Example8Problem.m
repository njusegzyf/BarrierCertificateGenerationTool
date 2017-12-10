function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example8Problem()

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [x1-x2;
     x1+x2];

eps = [0.1,0.1];

% Constructing the theta constraint
theta1 = 2*x1-5;
theta2 = -x2;
theta3 = x2;
g_theta = [theta1,theta2];

% Constructing the psy constraint
psy1 = x1/4;
psy2 = x2/4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1/2;
g_zeta = [zeta1];

end
