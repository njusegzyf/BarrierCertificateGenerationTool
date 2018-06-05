function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftTSC_C7Problem()

% independent variables
syms x1 x2;
vars = [x1, x2];

% Constructing the vector field dx/dt = f
f = [-x1 + x1*x2;
     -x2];

eps = [0.001, 0.001];

% Constructing the theta constraint
g_theta = vars + 2/3;

% Constructing the psy constraint
g_psy = (vars + 2)/4;

% Constructing the zeta constraint
zeta1 = x1 + 1/2;
zeta2 = x2 - 1/2;
g_zeta = [zeta1, zeta2];

end
