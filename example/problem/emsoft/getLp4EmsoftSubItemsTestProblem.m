function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftSubItemsTestProblem()

% independent variables
syms x1 x2 x3 x4 x5 x6 x7;
vars = [x1, x2, x3, x4, x5, x6, x7];

% Constructing the vector field dx/dt = f
f = ones(7, 1);

eps = [0, 0];

% Constructing the theta constraint
% g_theta = (vars+2)/4;
g_theta = vars/2;

% Constructing the psy constraint
g_psy = (vars+1)/2;

% Constructing the zeta constraint
g_zeta = (vars-1.8)*5;

end
