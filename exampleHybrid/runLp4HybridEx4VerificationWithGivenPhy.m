function [lpVer, solveResVer, resNorms] = runLp4HybridEx4VerificationWithGivenPhy()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx4Problem();
x1 = vars(1);
x2 = vars(2);
x3 = vars(3);

% Set the degree of lambda and re
pLambdaDegree = 2;
pReDegree = 2;

phys = [(6337886250869*x1^2)/17592186044416 + (16208047838125*x1*x2)/17592186044416 - (7537036577*x1*x3)/4398046511104 + (71222733821139*x1)/2199023255552 + (2269054183117*x2^2)/4398046511104 - (10488280259*x2*x3)/8796093022208 + (133629312460281*x2)/4398046511104 + (996198259*x3^2)/2199023255552 - (472501718879*x3)/17592186044416 + 2198450678280699/4398046511104,...
    (349482613*x2)/8796093022208 - (4286304659*x1)/17592186044416 + (1361162529*x3)/17592186044416 - (6983509*x1*x2)/8796093022208 - (6253025*x1*x3)/17592186044416 - (10133639*x2*x3)/17592186044416 - (2498639269*x1^2)/8796093022208 - (22944611*x2^2)/17592186044416 - (19114477*x3^2)/8796093022208 - 8854019971/2199023255552];



lpVer = lp4.HybridLinearProgramVerificationWithGivenPhy.createWithRou(vars, stateNum,...
    fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree,...
    phys);

[lpVer, solveResVer, resNorms] = lpVer.solveWithLinprog();
lp4.Lp4Config.displaySolveRes(solveResVer, resNorms);

warning('on')
echo off;

end
