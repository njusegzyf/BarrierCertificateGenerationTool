function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx7Problem()

% collin

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [x^2 + 2 * x * y + 3 * y^2;
     2 * y * (2 * x + y)];

eps = [0.0001, 0.0001];

% Constructing the theta constraint
% Note: use theta
% g_theta = [x + 2.5, y - 0.5] or g_theta = [(x + 2.4)/0.8, (y - 0.6)/0.8] 
g_theta = [(x + 2.4)/0.8, (y - 0.6)/0.8];

% Constructing the psy constraint
psy1 = (x + 3) / 6;
psy2 = (y + 3) / 6;
g_psy = [psy1, psy2];

% Constructing the zeta constraint
zeta1 = x - 1;
zeta2 = y - 1;
g_zeta = [zeta1, zeta2];

end
