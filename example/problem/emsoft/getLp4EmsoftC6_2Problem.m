function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftC6_2Problem()

%====
%B=69.46478566+94.0516*x1+83.3774*x2-20.3721*x1^2-6.6288*x2^2+32.0463*x1*x2;
%====

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [x2; -x1+1/3*x1^3-x2];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = x1-1;
theta2 = x2+0.5;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x1+1.5);
zeta2 = (x2+1.5);
g_zeta = [zeta1, zeta2];

end
