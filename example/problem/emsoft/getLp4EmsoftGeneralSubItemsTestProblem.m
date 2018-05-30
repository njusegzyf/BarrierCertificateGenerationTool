function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftGeneralSubItemsTestProblem(varNum)

% independent variables
vars = sym('x', [1, varNum]);

% Constructing the vector field dx/dt = f
f = ones(varNum, 1);

eps = [0, 0];

% Constructing the theta constraint
g_theta = (vars+0.1)*5;

% Constructing the psy constraint
g_psy = (vars+2)/4;

% Constructing the zeta constraint
g_zeta = (vars-1.8)*5;

end
