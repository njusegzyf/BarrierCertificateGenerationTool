function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Emsoft6d1Problem()

%- 0.1725*x1^2 - 0.6077*x1*x2 - 0.5707*x1*x3 - 0.03331*x1*x4 - 0.47*x1*x5 - 0.03331*x1*x6 + 0.3315*x1 + 0.209*x2^2 + 0.2005*x2*x3 - 0.03296*x2*x4 - 0.03367*x2*x5 - 0.03296*x2*x6 + 0.5343*x2 + 2.997*x3^2 - 0.2995*x3*x4 - 0.3037*x3*x5 - 0.2995*x3*x6 + 5.1*x3 + 1.427*x4^2 - 0.3047*x4*x5 - 0.3164*x4*x6 + 1.916*x4 + 1.989*x5^2 - 0.3047*x5*x6 + 2.912*x5 + 1.427*x6^2 + 1.916*x6 + 6.718

% independent variables
syms x1 x2 x3 x4 x5 x6;
vars = [x1; x2; x3; x4; x5; x6];

% Constructing the vector field dx/dt = f
f = [x1*x3; x1*x5;(x4-x3)*x3-2*x5^2;-(x4-x3)^2+(x6^2-x1^2);x2*x6+(x3-x4)*x5;2*x2*x5-x3*x6];

eps = [0.0001,0.0001];

% Constructing the theta constraint
theta1 = x1-1;
theta2 = x2-1;
theta3 = x3-1;
theta4 = x4-1;
theta5 = x5-1;
theta6 = x6-1;
g_theta = [theta1,theta2,theta3,theta4,theta5,theta6];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
psy3 = (x3+2)/4;
psy4 = (x4+2)/4;
psy5 = (x5+2)/4;
psy6 = (x6+2)/4;
g_psy = [psy1, psy2, psy3, psy4, psy5, psy6];

% Constructing the zeta constraint
zeta1 = 2*(x1+1.25);
zeta2 = 2*(x2+1.25);
zeta3 = 2*(x3+1.25);
zeta4 = 2*(x4+1.25);
zeta5 = 2*(x5+1.25);
zeta6 = 2*(x6+1.25);
g_zeta = [zeta1,zeta2,zeta3,zeta4,zeta5,zeta6];

end
