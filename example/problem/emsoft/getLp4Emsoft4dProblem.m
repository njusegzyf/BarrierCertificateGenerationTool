function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Emsoft4dProblem()

%====
%B= 23.35649653+11.9570*x1+5.6366*x2+38.1902*x3+17.9025*x4-6.8542*x3*x4+4.8428*x1*x4-9.4084*x1*x3-50.6959*x1^2-20.7429*x2^2-52.4027*x3^2-49.9831*x4^2-2.9252*x1*x2+0.8951*x2*x3+20.4207*x2*x4;
%====

% independent variables
syms x1 x2 x3 x4;
vars = [x1, x2, x3, x4];

% Constructing the vector field dx/dt = f
f = [-x1+x2^3-3*x3*x4; -x1-x2^3; x1*x4-x3; x1*x3-x4^3];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = 2*x1;
theta2 = 2*x2;
theta3 = 2*x3;
theta4 = 2*x4;

g_theta = [theta1, theta2, theta3, theta4];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
psy3 = (x3+2)/4;
psy4 = (x4+2)/4;
g_psy = [psy1, psy2, psy3, psy4];

% Constructing the zeta constraint
zeta1 = 4/3*(x1+1);
zeta2 = 4/3*(x2+1);
zeta3 = 4/3*(x3+1);
zeta4 = 4/3*(x4+1);
g_zeta = [zeta1, zeta2, zeta3, zeta4];

end
