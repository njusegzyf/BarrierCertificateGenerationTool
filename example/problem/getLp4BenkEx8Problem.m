function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx8Problem()

% dai

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [-2 * x + x^2 + y;
     x - 2 * y + y^2];

eps = [0.0001, 0.0001];

% Note: use
% phy = [(x + 1)/2, (y + 1)/2], zeta = [2*x - 1, 2*y - 1] or
% phy = [(x + 2)/4, (y + 2)/4], zeta = [2*x - 2, 2*y - 2]

% Constructing the theta constraint
theta1 = 5 * (x + 0.1);
theta2 = 5 * (y + 0.1);
g_theta = [theta1, theta2];

% Constructing the psy and zeta constraint
g_psy = [(x + 1)/2, (y + 1)/2];
g_zeta = [2*x - 1, 2*y - 1];

end
