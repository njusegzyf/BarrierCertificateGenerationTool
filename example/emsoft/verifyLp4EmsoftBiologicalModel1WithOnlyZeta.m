function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftBiologicalModel1WithOnlyZeta()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, ~, ~, zeta] = getLp4EmsoftBiologicalModel1();

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);
x5 = vars(5);
x6 = vars(6);
x7 = vars(7);

% initPhy = x1 + x2 + x3 + x4 + x5 + x6 + x7;
% this is computed from above initPhy with degree = 1
% initPhy = (3806265116065395*x1)/1152921504606846976 + (7069502577718175*x2)/2305843009213693952 + (3072282078226279*x3)/1152921504606846976 - (5857062640944425*x5)/9223372036854775808 + (5845620132802913*x6)/9223372036854775808 - 4962735518231069/576460752303423488;
 
initPhy = -2.549*x3-2.995*x6-3.576*x4-2.557*x1-3.389*x5-1.348*x2-2.615*x7+.754*x4^2+0.5552e-1*x4*x1+.7386*x5^2+0.4643e-1*x5*x7+.1847*x5*x2-.7219*x6*x3+.7259*x6^2+.3867*x1^2+1.299*x1*x2+.1235*x1*x3-.6037*x1*x5+.1977*x1*x6-.4632*x1*x7+0.4123e-1*x4*x2-.713*x4*x3+.7813*x4*x5+.215*x4*x6-.517*x4*x7-2.273*x2^2+1.446*x2*x3-.1724*x2*x6-0.3034e-1*x2*x7+0.9971e-1*x5*x3-.7006*x5*x6+0.9994e-1*x7*x3+.6089*x7*x6+.5814*x7^2+.4251*x3^2+22.71;

pLambdaDegree = 0;



% run and verify
lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, [], [], zeta, pLambdaDegree, initPhy);
[lpVer, solveResVer, resNorms] = lpVer.solve();

lp4.Lp4Config.printVerifyWithOnlyZetaResult(solveResVer, resNorms);

warning('on')

echo off;
end
