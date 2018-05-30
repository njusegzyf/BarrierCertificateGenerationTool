function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftSubItemsTestT7()

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

t1=-0.6613896068+0.0751*x1+0.0980*x2+0.0627*x3+0.0733*x4+0.0831*x5+0.1089*x6+0.0926*x7;
t2=0.01248137432-0.1466*x1+0.3543*x2+0.0496*x3-0.0484*x4+0.0021*x5-0.0837*x6-0.0656*x7;
t3=-0.03371294344+0.0679*x1-0.0486*x2+0.0688*x3+0.0646*x4+0.0732*x5-0.0640*x6-0.3324*x7;
t4=0.006759812806-0.1551*x1-0.0049*x2-0.1684*x3+0.0651*x4+0.1203*x5+0.1637*x6-0.0589*x7;
t5=-0.02408465683+0.0648*x1+0.0287*x2-0.0727*x3-0.1692*x4-0.1448*x5+0.1287*x6-0.0931*x7;
t6=-0.01954587036-0.0185*x1+0.0151*x2-0.0749*x3+0.1658*x4-0.1776*x5-0.0359*x6-0.0195*x7;
t7=0.02000802306-0.0817*x1-0.0253*x2+0.1753*x3+0.0379*x4-0.0591*x5+0.1266*x6-0.0088*x7;
t8=-0.05256722291-0.1226*x1-0.0781*x2-0.0054*x3-0.0696*x4-0.0339*x5-0.0750*x6-0.0160*x7;

testItems = [t1, t2, t3, t4, t5, t6, t7, t8];



% run test on each t
for i = 1 : 8
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Test Item t7', int2str(i)]);
    lp4.Lp4Config.displayDelimiterLine();
    
    % run and verify
    lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, [], [], pLambdaDegree, testItems(i)^2);
    [lpVer, solveResVer, resNorms] = lpVer.solve();

    lp4.Lp4Config.printVerifyWithOnlyThetaResult(solveResVer, resNorms);
end

warning('on')

echo off;
end
