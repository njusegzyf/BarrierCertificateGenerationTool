function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftSingleVarTest2Problem()

% independent variables
syms x;
vars = [x];

% Constructing the vector field dx/dt = f
f = [1];

eps = [0, 0];

% Constructing the theta constraint
g_theta = [x];

% Constructing the psy constraint
g_psy = (vars+1)/2;

% Constructing the zeta constraint
g_zeta = (vars-1.8)*5;

end
