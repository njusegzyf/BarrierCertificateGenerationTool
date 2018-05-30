function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftSubItemsTestT3()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, ~, ~] = getLp4EmsoftSubItemsTestProblem();

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);
x5 = vars(5);
x6 = vars(6);
x7 = vars(7);

pLambdaDegree = 0;

t1=-0.7979222438+0.0711*x1+0.0738*x2+0.0333*x3+0.2435*x4+0.0662*x5+0.0512*x6+0.0652*x7;
t2=0.002896807087-0.1160*x1+0.3373*x2+0.1622*x3-0.0620*x4+7.1097e-04*x5-0.0685*x6-0.0180*x7;
t3=-0.007751260181-0.1088*x1+0.1280*x2-0.3236*x3-0.0113*x4+0.0099*x5+0.0482*x6+0.0387*x7;
t4=0.02711279985-0.0584*x1-0.0131*x2-0.0276*x3+0.1587*x4+0.0948*x5-0.1299*x6-0.1624*x7;
t5=0.01941396956-0.1578*x1-0.0450*x2+0.0639*x3+0.0418*x4+0.1130*x5+0.1510*x6+0.0388*x7;
t6=-0.01887303744+0.0843*x1+0.0170*x2-0.0165*x3-0.1200*x4+0.1653*x5+0.0324*x6-0.0788*x7;
t7=0.0119447449+0.0397*x1+0.0391*x2+0.0025*x3+0.0411*x4-0.0742*x5+0.1253*x6-0.1193*x7;
t8=-0.03377491272-0.0579*x1-0.0354*x2-3.7089e-04*x3-0.0584*x4-0.0332*x5-0.0133*x6-0.0477*x7;

testItems = [t1, t2, t3, t4, t5, t6, t7, t8];



% run test on each t
for i = 1 : 8
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Test Item t3', int2str(i)]);
    lp4.Lp4Config.displayDelimiterLine();
    
    % run and verify
    lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, [], [], pLambdaDegree, testItems(i)^2);
    [lpVer, solveResVer, resNorms] = lpVer.solve();

    lp4.Lp4Config.printVerifyWithOnlyThetaResult(solveResVer, resNorms);
end

warning('on')

echo off;
end
