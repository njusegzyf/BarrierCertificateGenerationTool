function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example8Problem()

% From FM2016 ref 05 Computation of parametric barrier functions for dynamical systems using interval analysis
% Example 6: (Lorenz) [25] We consider the Lorenz system with a limit cycle. 
% The initial region is taken inside the limit cycle and the unsafe outside of it.

% independent variables
syms x1 x2 x3;
vars = [x1 x2 x3];

% Constructing the vector field dx/dt = f
f = [10*(x2-x1); x1*(28-x3)-x2; x1*x2-8/3*x3];

eps = [0.1,0.1];

% Constructing the theta constraint
theta1 = x1+15;
theta2 = x2+15;
theta3 = (x3-10/5);
g_theta = [theta1,theta2,theta3];

% Constructing the psy constraint
psy1 = (x1+20)/40;
psy2 = (x2+20)/10;
psy3 = (x3+20)/40;
g_psy = [psy1, psy2, psy3];

% Constructing the zeta constraint
zeta1 = x1+17;
zeta2 = x2+15;
zeta3 = x3/5;

g_zeta = [zeta1,zeta2,zeta3];

end
