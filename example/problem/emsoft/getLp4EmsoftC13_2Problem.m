function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftC13_2Problem()

%====
%B=26.50266211-0.1484*x1-112.0672*x2+0.0508*x1*x2+45.2351*x1^2+5.4202*x2^2;
%====

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [-x1+x1*x2; -x2];

eps = [1, 1];

% Constructing the theta constraint
theta1 = x1+1.5;
theta2 = x2+1.5;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x1+0.5);
zeta2 = (x2-0.5);
g_zeta = [zeta1, zeta2];

end
