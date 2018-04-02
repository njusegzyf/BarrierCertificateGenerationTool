function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Emsoft4d2Problem()

%====
%B=5.782725478+34.8916*x1-72.8435*x2-3.6512*x3+70.1655*x4+5.1731*x1^2-4.3645*x2^2+2.3811*x3^2+4.6440*x2*x3-3.4586*x4^2-8.9548*x1*x2-2.9863*x1*x3+8.6539*x1*x4+8.6151*x2*x4-9.5623*x3*x4;
%====

% independent variables
syms x1 x2 x3 x4;
vars = [x1, x2, x3, x4];

% Constructing the vector field dx/dt = f
f = [-1/2*x1^2-1/2*x2^2-1/8*x3^2-2*x2*x3+2*x4^2+1; -x1*x2-1; -x1*x3; -x1*x4];

eps = [0.001, 0.001];

% Constructing the theta constraint
theta1 = 2*(x1-1.5);
theta2 = 2*(x2-1.5);
theta3 = 2*(x3-1.5);
theta4 = 2*(x4-1.5);

g_theta = [theta1, theta2, theta3, theta4];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
psy3 = (x3+2)/4;
psy4 = (x4+2)/4;
g_psy = [psy1, psy2, psy3, psy4];

% Constructing the zeta constraint
zeta1 = 2*(x1+2);
zeta2 = 2*(x2+2);
zeta3 = 2*(x3+2);
zeta4 = 2*(x4+2);
g_zeta = [zeta1, zeta2, zeta3, zeta4];

end
