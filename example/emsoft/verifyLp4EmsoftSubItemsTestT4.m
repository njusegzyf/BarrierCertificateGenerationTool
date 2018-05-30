function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftSubItemsTestT4()

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

t1=-0.6536865448+0.0589*x1+0.0942*x2+0.1415*x3+0.0452*x4+0.0938*x5+0.0781*x6+0.0858*x7;
t2=-0.02051429041+0.1240*x1-0.3365*x2-0.0308*x3+0.1543*x4+0.0207*x5+0.0711*x6+0.0102*x7;
t3=-0.02383259693+0.0965*x1-0.1035*x2+0.0052*x3-0.3163*x4-0.0499*x5+0.0308*x6+0.0504*x7;
t4=-0.001247182893+0.1035*x1-0.0210*x2+0.1616*x3+0.0084*x4-0.0295*x5-0.2044*x6-0.1101*x7;
t5=0.002657619175-0.0447*x1-0.0306*x2-0.0047*x3-0.0590*x4+0.2233*x5+0.0141*x6-0.1337*x7;
t6=-0.01264666133-0.0722*x1-0.0477*x2-0.0345*x3-0.0086*x4+0.0595*x5-0.1352*x6+0.1250*x7;
t7=-0.01270577416+0.1256*x1+0.0547*x2-0.1363*x3+0.0083*x4+0.0289*x5-0.0479*x6-0.0105*x7;
t8=-0.04555303796-0.0605*x1-0.0308*x2-0.0616*x3-0.0106*x4-0.0609*x5-0.0291*x6-0.0715*x7;

testItems = [t1, t2, t3, t4, t5, t6, t7, t8];



% run test on each t
for i = 1 : 8
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Test Item t4', int2str(i)]);
    lp4.Lp4Config.displayDelimiterLine();
    
    % run and verify
    lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, [], [], pLambdaDegree, testItems(i)^2);
    [lpVer, solveResVer, resNorms] = lpVer.solve();

    lp4.Lp4Config.printVerifyWithOnlyThetaResult(solveResVer, resNorms);
end

warning('on')

echo off;
end
