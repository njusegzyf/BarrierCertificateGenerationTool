function [lpVer, solveResVer, resNorms] = runLp4HybridEx44VerificationWithGivenPhy()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx44Problem();
x1 = vars(1);
x2 = vars(2);
x3 = vars(3);

% Set the degree of lambda and re
pLambdaDegree = 2;
pReDegree = 2;

phys = [(2965460925185*x1^2)/8796093022208 + (8431731108817*x1*x2)/8796093022208 - (205715891*x1*x3)/8796093022208 + (204620066094885*x1)/17592186044416 + (5907589313043*x2^2)/2199023255552 - (11619599963*x2*x3)/17592186044416 + (533220417383809*x2)/8796093022208 + (47102617*x3^2)/8796093022208 - (131667725313*x3)/17592186044416 + 4396825861060565/8796093022208,...
    (139670267*x2)/8796093022208 - (4830874815*x1)/17592186044416 + (339826119*x3)/8796093022208 - (6924753*x1*x2)/8796093022208 - (1696773*x1*x3)/8796093022208 - (2977751*x2*x3)/8796093022208 - (4764760887*x1^2)/17592186044416 - (9236145*x2^2)/17592186044416 - (19547611*x3^2)/17592186044416 - 46800231393/17592186044416];



lpVer = lp4.HybridLinearProgramVerificationWithGivenPhy.createWithRou(vars, stateNum,...
    fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree,...
    phys);

[lpVer, solveResVer, resNorms] = lpVer.solve();
lp4.Lp4Config.displaySolveRes(solveResVer, resNorms);

warning('on')
echo off;

end
