function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftSubItemsTestT5()

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

t1=-0.6707850144+0.0563*x1+0.0857*x2+0.0617*x3+0.0867*x4+0.0496*x5+0.1716*x6+0.0838*x7;
t2=0.0109963588-0.1413*x1+0.3559*x2+0.0417*x3-0.0497*x4+0.0129*x5-0.0739*x6-0.0167*x7;
t3=0.01446619004-0.0333*x1-0.0270*x2+0.0070*x3+0.0524*x4+0.3527*x5-0.0205*x6-0.0605*x7;
t4=-0.01995415779+0.1119*x1-0.0081*x2+0.2031*x3+0.0520*x4-0.0248*x5-0.1617*x6-0.0839*x7;
t5=0.003989698338-0.1224*x1-0.0171*x2-0.0154*x3+0.2167*x4-0.0533*x5-3.7420e-04*x6-0.0490*x7;
t6=0.004762724387-0.0610*x1-0.0345*x2+0.0563*x3+0.0111*x4+0.0158*x5-0.0623*x6+0.1795*x7;
t7=0.02083036875-0.0805*x1-0.0323*x2+0.1359*x3-0.0474*x4-0.0054*x5+0.1167*x6-0.0329*x7;
t8=0.04412029044+0.0714*x1+0.0474*x2+0.0224*x3+0.0573*x4+0.0107*x5+0.0643*x6+0.0430*x7;

testItems = [t1, t2, t3, t4, t5, t6, t7, t8];



% run test on each t
for i = 1 : 8
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Test Item t5', int2str(i)]);
    lp4.Lp4Config.displayDelimiterLine();
    
    % run and verify
    lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, [], [], pLambdaDegree, testItems(i)^2);
    [lpVer, solveResVer, resNorms] = lpVer.solve();

    lp4.Lp4Config.printVerifyWithOnlyThetaResult(solveResVer, resNorms);
end

warning('on')

echo off;
end
