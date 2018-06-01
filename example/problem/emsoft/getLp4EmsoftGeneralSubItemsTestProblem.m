function [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4EmsoftGeneralSubItemsTestProblem(varNum)
% To verify an expr can be expanded with bases, get this problem, create a LinearProgram4Verification3,
% set the expr to be the phy and the bases to be the theta, and then verify the problem. 

% independent variables
vars = sym('x', [1, varNum]);

% Constructing the vector field dx/dt = f
f = ones(varNum, 1);

eps = [0, 0];

% Constructing the theta constraint
g_theta = [];

% Constructing the psy constraint
g_psy = [];

% Constructing the zeta constraint
g_zeta = [];

end
