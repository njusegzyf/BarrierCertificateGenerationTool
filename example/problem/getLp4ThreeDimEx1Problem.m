function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4ThreeDimEx1Problem()

syms x1 x2 x3;
vars = [x1, x2, x3];

f = [-x2;
    -x3;
    -x1 - 2 * x2 - x3 + x1^3];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
theta1 = x1+0.25;
theta2 = x2+0.25;
theta3 = x3+0.75;
g_theta = [theta1,theta2,theta3];

% Constructing the psy constraint
psy1 = (x1+2)/4;
psy2 = (x2+2)/4;
psy3 = (x3+2)/4;
g_psy = [psy1, psy2, psy3];

% Constructing the zeta constraint
zeta1 = x1-1;
zeta2 = x2+2;
zeta3 = x3+2;
g_zeta = [zeta1,zeta2,zeta3];

end
