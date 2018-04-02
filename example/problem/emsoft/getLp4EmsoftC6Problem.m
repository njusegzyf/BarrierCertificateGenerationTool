function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftC6Problem()

%====
% Init=(x1-1.5)^2+x2^2-0.25;
% Unsafe=(x1+1)^2+(x2+1)^2-0.64;
% Inv1=(x1+2)*(x1-2);
% Inv2=(x2-2)*(x2+2);
%B=61.94063789+90.6204*x1+81.3145*x2-22.6008*x1^2-7.0696*x2^2+31.0851*x1*x2;
%====

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [x2; -x1+1/3*x1^3-x2];

eps = [0.0001, 0.0001];

% 原始约束
% % Constructing the theta constraint
% theta1 = 4*((x1-1.5)^2+x2^2-0.25);
% g_theta = [theta1];
% 
% % Constructing the psy constraint
% psy1 = (x1+2)/4;
% psy2 = (x2+2)/4;
% g_psy = [psy1, psy2];
% 
% % Constructing the zeta constraint
% zeta1 = 25/16*((x1+1)^2+(x2+1)^2-0.64);
% g_zeta = [zeta1];

% 修改后的约束
% Constructing the theta constraint
theta1 = 4*((x1-1.5)^2+x2^2-0.25);
theta2 = x1-1;
theta3 = x2+0.5;
g_theta = [theta1, theta2, theta3];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = 25/16*((x1+1)^2+(x2+1)^2-0.64);
zeta2 = (x1+1.8)/1.6;
zeta3 = (x2+1.8)/1.6;
g_zeta = [zeta1, zeta2, zeta3];

end
