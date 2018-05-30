function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftSubItemsTestT6()

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

t1=-0.7468276997+0.0473*x1+0.0740*x2+0.0419*x3+0.0636*x4+0.2165*x5+0.0446*x6+0.0979*x7;
t2=0.006142059864-0.1096*x1+0.3190*x2+0.0578*x3-0.0398*x4-0.0062*x5-0.2033*x6-0.0342*x7;
t3=-0.01901876442+0.1048*x1-0.1610*x2+0.0419*x3+0.0332*x4+0.0333*x5-0.2965*x6-0.0520*x7;
t4=0.03242018425-0.1181*x1-0.0207*x2-0.0609*x3+0.0942*x4+0.1843*x5+0.0125*x6-0.1282*x7;
t5=-0.003518515156-0.1012*x1-0.0360*x2-0.1648*x3-0.0310*x4+0.0068*x5-0.0694*x6+0.1566*x7;
t6=0.007787406776+0.0621*x1+0.0075*x2-0.0435*x3-0.1906*x4+0.0925*x5+0.0077*x6-0.0420*x7;
t7=-0.015088187+0.0899*x1+0.0523*x2-0.1365*x3+0.0518*x4-0.0491*x5-0.0036*x6-0.0631*x7;
t8=0.03312405787+0.0529*x1+0.0286*x2+0.0086*x3+0.0327*x4+0.0560*x5+0.0024*x6+0.0557*x7;

testItems = [t1, t2, t3, t4, t5, t6, t7, t8];



% run test on each t
for i = 1 : 8
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Test Item t6', int2str(i)]);
    lp4.Lp4Config.displayDelimiterLine();
    
    % run and verify
    lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, [], [], pLambdaDegree, testItems(i)^2);
    [lpVer, solveResVer, resNorms] = lpVer.solve();

    lp4.Lp4Config.printVerifyWithOnlyThetaResult(solveResVer, resNorms);
end

warning('on')

echo off;
end
