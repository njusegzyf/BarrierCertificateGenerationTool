function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftSubItemsTestT1()

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

t11=-0.6477873095+0.0645*x1+0.0897*x2+0.0950*x3+0.0849*x4+0.0772*x5+0.0790*x6+0.1042*x7;
t12=0.006865717149-0.2835*x1+0.3063*x2+0.0090*x3-0.0335*x4+0.0153*x5-0.0420*x6-0.0058*x7;
t13=-0.01499333122-0.2001*x1-0.1534*x2-0.1403*x3+0.0696*x4+0.0557*x5+0.1540*x6+0.0759*x7;
t14=0.001667727555+0.0673*x1+0.0721*x2-0.1051*x3-0.1645*x4-0.1315*x5+0.1122*x6+0.1489*x7;
t15=-0.01965350658-0.0975*x1-0.0931*x2+0.1402*x3-0.0056*x4-0.1766*x5-0.0640*x6+0.0745*x7;
t16=-0.01604063189-0.0401*x1-0.0638*x2+0.0803*x3-0.1955*x4+0.1083*x5+0.0193*x6-0.0290*x7;
t17=-0.02248260921-0.0048*x1-0.0306*x2-0.0979*x3-0.0304*x4+0.0553*x5-0.1853*x6+0.1032*x7;
t18=0.06039745294+0.0252*x1+0.0272*x2+0.0895*x3+0.0452*x4+0.0800*x5+0.0317*x6+0.1347*x7;

testItems = [t11, t12, t13, t14, t15, t16, t17, t18];



% run test on each t
for i = 1 : 8
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Test Item t1', int2str(i)]);
    lp4.Lp4Config.displayDelimiterLine();
    
    % run and verify
    lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, [], [], pLambdaDegree, testItems(i)^2);
    [lpVer, solveResVer, resNorms] = lpVer.solve();

    lp4.Lp4Config.printVerifyWithOnlyThetaResult(solveResVer, resNorms);
end

warning('on')

echo off;
end
