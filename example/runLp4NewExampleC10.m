function [lpVer, solveResVer, resNorms] = runLp4NewExampleC10() 

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, theta, psy, zeta] = getLp4NewExampleC10Problem();

% Set the degree of phy
degree = 4;

lambda = 0;



% solve the problem
lpVer = lp4.LinearProgram4Verification2.createWithRou(vars, f, eps, theta, psy, zeta, degree, lambda, 0);
[lpVer, solveResVer, resNorms] = lpVer.solve();
    
warning('on')

echo off;
end
