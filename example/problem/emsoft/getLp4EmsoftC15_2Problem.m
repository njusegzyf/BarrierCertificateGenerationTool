function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftC15_2Problem()

%====
%B=-44.80572631+44.6479*x1+71.3698*x2+25.8684*x1*x2-11.3483*x1^2-18.2959*x2^2;
%====

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [x2+2*x1*x2; -x1+2*x1^2-x2^2];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = x1;
theta2 = x2-1;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x1+1);
zeta2 = (x2+0.5);
g_zeta = [zeta1, zeta2];

end
