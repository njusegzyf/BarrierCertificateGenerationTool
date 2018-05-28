function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftTest1Problem()

% independent variables
syms x y;
vars = [x, y];

% Constructing the vector field dx/dt = f
f = [-x;
    -y];

eps = [0.001, 0.001];

% Constructing the theta constraint
g_theta = (vars+0.1)*5;

% Constructing the psy constraint
g_psy = (vars+2)/4;

% Constructing the zeta constraint
g_zeta = (vars-1.8)*5;

end
