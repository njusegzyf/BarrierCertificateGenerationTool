function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example9Problem()

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [-x1+x1*x2;
     x2-x1*x2];

eps = [0.1,0.1];

% Constructing the theta constraint
theta1 = 625/4*((x1-0.8)^2+(x2-0.2)^2);
theta2 = x1-0.72;
theta3 = x2-0.12;
g_theta = [theta1,theta2,theta3];

% Constructing the psy constraint
psy1 = 5/4*x1-1/8;
psy2 = 5/4*x2-1/8;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 100*(x1-0.6)^2+100*(x2-0.6)^2;
zeta2 = x1-0.5;
zeta3 = x2-0.5;
g_zeta = [zeta1,zeta2,zeta3];

end
