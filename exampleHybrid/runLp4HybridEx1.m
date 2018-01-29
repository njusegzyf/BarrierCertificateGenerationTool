function [hlp, solveRes] = runLp4HybridEx1() 

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards] = getLp4HybridEx1Problem();


% Set the degree of phy, lambda and re
degree = 2;
pLambdaDegree = 2;
pReDegree = 2;


import lp4util.Partition
phyRange = Partition(-1, 1);
pLambdaRange = Partition(-1, 1);
pReRange = Partition(-1, 1);
phyRangeInVerify = Partition(-10, 10);



import lp4.HybridLinearProgram
hlp = HybridLinearProgram.create(vars, stateNum,...
fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta,guards, degree, pLambdaDegree, pReDegree,...
phyRange, pLambdaRange, pReRange);

[hlp, solveRes] = hlp.solve();

% warning('on')

echo off;
end
