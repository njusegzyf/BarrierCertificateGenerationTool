function res = runLp4EmsoftC5WithFmincom() 

% C5

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

[vars, f, ~, g_theta, g_psy, g_zeta] = getLp4BenkEx61Problem();

eps = [0.1, 0.1];

% Set the degree of phy
degree = 2;
pLambdaDegree = 2;



res = lp4.solveBlpDirectlyWithFmincon(vars, f, eps, g_theta, g_psy, g_zeta, degree, pLambdaDegree);

warning('on')

echo off;
end
