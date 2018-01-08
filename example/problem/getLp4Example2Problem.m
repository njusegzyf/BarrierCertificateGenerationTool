function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example2Problem()

% From FM2016 ref 16 A Framework for Worst-Case and Stochastic Safety Verification Using Barrier Certificates
% V. EXAMPLES A. Continuous System

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [x2;
    -x1+(1/3)*x1^3-x2];

eps = [0.00001,0.00001];

% Constructing the theta constraint
%theta1 = 4*(x1-1.5)^2+4*x2^2;
theta1 = x2+1/2;
theta2 = x1-1;
g_theta = [theta1,theta2];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
%zeta1 = 6.25*(x1+1)^2+6.25*(x2+1)^2;
zeta1 = (x1+1.3)/0.6;
zeta2 = (x2+1.3)/0.6;
g_zeta = [zeta1,zeta2];

end
