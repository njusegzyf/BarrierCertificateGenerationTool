function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example3Problem()

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [2*x1-x1*x2;
     2*x1^2-x2;];

eps = [0.001,0.001];

% Constructing the theta constraint
theta1 = x1^2+(x2+2)^2;

g_theta = [theta1];

% Constructing the psy constraint
psy1 = (x1+7)/14;
psy2 = (x2+7)/14;

g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x1^2+(x2-5.2)^2;

g_zeta = [zeta1];

end
