function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftSubItemsTestT2()

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

t1=-0.4328085014-0.0094*x1+0.2331*x2+0.0718*x3+0.0395*x4+0.0617*x5+0.0478*x6+0.0686*x7;
t2=0.1662614695-0.1385*x1+0.3577*x2+0.0031*x3-0.0769*x4-0.0388*x5-0.0844*x6-0.0510*x7;
t3=-0.00928146193+0.0724*x1-0.0085*x2+0.1485*x3+0.0472*x4+0.0169*x5-0.1502*x6-0.1129*x7;
t4=0.0114258275-0.0999*x1-0.0096*x2-0.0432*x3+0.1066*x4+0.1599*x5+0.0262*x6-0.0873*x7;
t5=-0.003067590861+0.0692*x1+0.0088*x2-0.0094*x3-0.1383*x4+0.0944*x5+0.0366*x6-0.0608*x7;
t6=0.01249758249+0.0488*x1+0.0310*x2+0.0180*x3+0.0502*x4-0.0645*x5+0.1143*x6-0.0891*x7;
t7=0.01704080872-0.0773*x1-0.0276*x2+0.1280*x3-0.0366*x4+0.0116*x5+0.0607*x6+0.0252*x7;
t8=-0.0368949996-0.0533*x1-0.0272*x2-0.0207*x3-0.0348*x4-0.0430*x5-0.0156*x6-0.0565*x7;

testItems = [t1, t2, t3, t4, t5, t6, t7, t8];



% run test on each t
for i = 1 : 8
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Test Item t2', int2str(i)]);
    lp4.Lp4Config.displayDelimiterLine();
    
    % run and verify
    lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, [], [], pLambdaDegree, testItems(i)^2);
    [lpVer, solveResVer, resNorms] = lpVer.solve();

    lp4.Lp4Config.printVerifyWithOnlyThetaResult(solveResVer, resNorms);
end

warning('on')

echo off;
end
