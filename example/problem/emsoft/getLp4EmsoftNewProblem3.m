function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftNewProblem3()

% Exp10

% independent variables
syms x y z v w;
vars = [x, y, z, v, w];

% Constructing the vector field dx/dt = f
f = [-0.1*x^2 - 0.4*x*v - x + y + 3*z + 0.5*v;
    y^2 - 0.5*y*w + x + z;
    0.5*z^2 + x - y + 2*z + 0.1*v - 0.5*w;
    -0.5*y + 0.5*z - 1.4*v - 1.7*w - 1.5*x; % y + 2*z + 0.1*v - 0.2*w - 1.5*x - 1.5*y - 1.5*z - 1.5*v - 1.5*w
    -0.5*z - 1.6*v - 1.5*x - 1.5*y - 1.5*w]; % z - 0.1*v - 1.5*x - 1.5*y - 1.5*z - 1.5*v - 1.5*w

eps = [0.001, 0.001];

% Constructing the theta constraint
% theta1 = x;
% theta2 = y;
% theta3 = z;
% theta4 = v;
% theta5 = w;
% g_theta = [theta1, theta2, theta3, theta4, theta5];
g_theta = vars;

% Constructing the psy constraint
% psy1 = (x+2)/4;
% psy2 = (y+2)/4;
% psy3 = (z+2)/4;
% psy4 = (v+2)/4;
% psy5 = (w+2)/4;
% g_psy = [psy1, psy2, psy3, psy4, psy5];
g_psy = (vars+2)/4; 

% Constructing the zeta constraint
% zeta1 = 2*(x+2);
% zeta2 = 2*(y+2);
% zeta3 = 2*(z+2);
% zeta4 = 2*(v+2);
% zeta5 = 2*(w+2);
% g_zeta = [zeta1, zeta2, zeta3, zeta4, zeta5];
g_zeta = 2*(vars+2);

end
