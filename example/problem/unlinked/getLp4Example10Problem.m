function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example10Problem()

% it has same f with example 3 but different constraints

% independent variables
syms x y;
vars = [x y];

% Constructing the vector field dx/dt = f
f = [1+x^2*y-2.5*x;
     1.5*x-x^2*y];

eps = [0.01, 0.01];

% Constructing the theta constraint
theta1 = 10*x-9;
theta2 = 10*y;
g_theta = [theta1,theta2];

% Constructing the psy constraint
psy1 = x/2;
psy2 = y/2;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 5*(x-1.6);
zeta2 = 5*(y-1.2);
g_zeta = [zeta1,zeta2];

end
