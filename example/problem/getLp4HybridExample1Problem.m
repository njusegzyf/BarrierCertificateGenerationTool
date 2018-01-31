function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4HybridExample1Problem()
% 该例子是将混成系统的例子转为连续系统的例子

% independent variables
syms x1 x2 x3;
vars = [x1 x2 x3];

% Constructing the vector field dx/dt = f
f = [x2;
    -x1+x3;
    -x1-2*x2-3*x3];

eps = [0.001,0.001];

% Constructing the theta constraint
theta1 = (x1^2+0.01*x2^2+0.01*x3^2-0.09)/0.02;
theta2 = (x1+1.01)/(2.02);
theta3 = (0.1*x2+1.01)/(2.02);
theta4 = (0.1*x3+1.01)/(2.02);
g_theta = [theta1, theta2, theta3, theta4];

% Constructing the psy constraint
psy1 = (x1+5.1)/10.2;
psy2 = (x2+20)/(40);
psy3 = (x3+20)/(40);
psy4 = (x1^2+x2^2+x3^2-0.03)/400;
g_psy = [psy1, psy2, psy3, psy4];

% Constructing the zeta constraint
zeta1 = (x2+20)/(40);
zeta2 = (x3+20)/(40);
zeta3 = (x1-5)*10;
g_zeta = [zeta1, zeta2, zeta3];

end
