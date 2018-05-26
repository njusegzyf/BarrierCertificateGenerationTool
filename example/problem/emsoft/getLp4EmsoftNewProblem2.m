function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftNewProblem2()

% Exp10

% independent variables
syms x y z v;
vars = [x, y, z, v];

% Constructing the vector field dx/dt = f
f = [y;
    -9.8*z + 1.6*z^3 + x*v^2;
    v;
    -6*x];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = x;
theta2 = y;
theta3 = z;
theta4 = v;

g_theta = [theta1, theta2, theta3, theta4];

% Constructing the psy constraint
psy1 = (x+2)/4;
psy2 = (y+2)/4;
psy3 = (z+2)/4;
psy4 = (v+2)/4;
g_psy = [psy1, psy2, psy3, psy4];

% Constructing the zeta constraint
zeta1 = 2*(x+2);
zeta2 = 2*(y+2);
zeta3 = 2*(z+2);
zeta4 = 2*(v+2);
g_zeta = [zeta1, zeta2, zeta3, zeta4];

end
