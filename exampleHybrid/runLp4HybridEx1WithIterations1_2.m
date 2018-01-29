function [lpVer, solveResVer, resNorms] = runLp4HybridEx1WithIterations1_2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx1Problem();
T = vars(1);
t = vars(2);

% Set the degree
degree = 2;
pLambdaDegree = 1;
pReDegree = 1;

maxIterations = 20;

initPhys = [(1805148783081117*t^2)/36893488147419103232 + 6123377283986903/2251799813685248,...
    (1733886359069771*T)/72057594037927936 + (4473506490099937*t)/144115188075855872 - (1413697157682783*T^2)/20769187434139310514121985316880384 - (8355650283930911*t^2)/73786976294838206464 - 4848994956863559/72057594037927936,...
    (8645256203897683*T^2)/83076749736557242056487941267521536 - (5781680151835247*T)/73786976294838206464 - 7039380949301427/72057594037927936];



import lp4.runAndVerifyHLPWithIterations1
[lpVer, solveResVer, resNorms] = runAndVerifyHLPWithIterations1(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initPhys);

warning('on')
echo off;

end
