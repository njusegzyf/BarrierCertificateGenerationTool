function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftC13Problem()

%====
%Init=(x1+1)^2+(x2+1)^2-0.25;
%Unsafe=(x1)^2+(x2-1)^2-0.25;
% Inv1=(x1+2)*(x1-2);
% Inv2=(x2-2)*(x2+2);
%B=20.32969076-42.3308*x1-118.3744*x2+19.8331*x1*x2+1.9158*x1^2-2.8513*x2^2;
%====

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [-x1+x1*x2; -x2];

eps = [1, 1];

% Constructing the theta constraint
theta1 = 4*((x1+1)^2+(x2+1)^2-0.25);
theta2 = x1+1.5;
theta3 = x2+1.5;
g_theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 4*((x1)^2+(x2-1)^2-0.25);
zeta2 = (x1+0.5);
zeta3 = (x2-0.5);
g_zeta = [zeta1, zeta2, zeta3];

end
