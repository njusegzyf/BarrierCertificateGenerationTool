function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example5Problem()

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [x2;
    2*x1-x2-x1^2*x2-x1^3];

eps = [0.1,0.1];

% Constructing the theta constraint
theta1 = 25/4*((x1+1)^2+(x2-2)^2);
theta2 = x1+7/5;
theta3 = x2-8/5;
g_theta = [theta1,theta2,theta3];

% Constructing the psy constraint
psy1 = (x1+3)/6;
psy2 = (x2+3)/6;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 25*((x1-1)^2+x2^2);
zeta2 = x1-4/5;
zeta3 = x2+1/5;
g_zeta = [zeta1,zeta2,zeta3];

end
