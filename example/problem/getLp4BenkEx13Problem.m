function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx13Problem()

% darbowx

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [x*(1-x^2-y^2)+y*((-1+x^2)^2+y^2);
     y*(1-x^2-y^2)-y*((-1+x^2)^2+y^2)];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = 6*x + 3;
theta2 = 2*y + 1;
g_theta = [theta1, theta2];

% Constructing the psy constraint
psy1 = (x + 1) / 2;
psy2 = (y + 1) / 2;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x;
zeta2 = y;
g_zeta = [zeta1, zeta2];

end
