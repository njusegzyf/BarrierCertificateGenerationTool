function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example6Problem()

% From FM2016 ref 24 Stabilization of polynomial dynamical systems using linear programming based on Bernstein polynomials
% Example 6

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

% independent variables
syms x1 x2 x3;
vars = [x1 x2 x3];

% Constructing the vector field dx/dt = f
f = [-x1+x2-x3; -x1*(x3+1)-x2; -x1+1.77*x1-4.8*x3];

eps = [0.001,0.001];

% Constructing the theta constraint
theta1 = x1-1;
theta2 = x2+0.5;
theta3 = x3+0.5;
g_theta = [theta1,theta2,theta3];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
psy3 = (x3+2)/4;
g_psy = [psy1, psy2, psy3];

% Constructing the zeta constraint
zeta1 = x1+1.5;
zeta2 = x2+1.5;
zeta3 = x3-0.5;

g_zeta = [zeta1,zeta2,zeta3];

end
