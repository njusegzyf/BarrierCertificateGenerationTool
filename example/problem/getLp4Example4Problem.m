function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example4Problem()

% independent variables
syms x1 x2;
vars = [x1 x2];

% Constructing the vector field dx/dt = f
f = [x1^2-2*x1+x2;
     x1+x2^2-2*x2;];

eps = [0.1,0.1];

% Constructing the theta constraint
theta1 = 100*(x1^2+x2^2);
g_theta = [theta1];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = (x1^2+x2^2-0.25)/7.75;
g_zeta = [zeta1];

end
