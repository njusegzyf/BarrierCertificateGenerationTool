function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftBiologicalModel2WithOnlyPsy()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, ~, psy, ~] = getLp4EmsoftBiologicalModel2();

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);
x5 = vars(5);
x6 = vars(6);
x7 = vars(7);
x8 = vars(8);
x9 = vars(9);

% initPhy = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9;
% this is computed from above initPhy with degree = 1
% initPhy = (1446247084809473*x2)/4611686018427387904 + (7231374831680253*x3)/2361183241434822606848 + (940704830558181*x4)/2305843009213693952 + (8317549420363109*x5)/4722366482869645213696 + (1807843820267583*x6)/590295810358705651712 + (8866204539489127*x7)/36893488147419103232 + (2834607139113733*x8)/18446744073709551616 + (8344624299873613*x9)/73786976294838206464 - 2547804529151789/2305843009213693952;

initPhy = 0.4788*x1^2 + 0.5227*x1*x2 + 2.492*x1*x3 + 0.3089*x1*x4 + 1.391*x1*x5 - 0.114*x1*x6 - 0.8628*x1*x7 + 0.2594*x1*x8 + 0.1612*x1*x9 - 8.967*x1 + 0.5521*x2^2 + 0.8764*x2*x3 + 0.9352*x2*x4 - 0.3628*x2*x5 - 0.04419*x2*x6 + 0.4258*x2*x7 - 0.08721*x2*x8 - 0.004147*x2*x9 - 6.376*x2 - 0.1785*x3^2 + 1.022*x3*x4 - 0.05612*x3*x5 - 0.002406*x3*x6 - 1.058*x3*x7 + 0.197*x3*x8 + 0.2396*x3*x9 - 10.73*x3 + 0.2405*x4^2 - 0.1054*x4*x5 + 0.3233*x4*x6 + 1.207*x4*x7 - 0.4179*x4*x8 - 0.115*x4*x9 - 6.947*x4 - 0.4099*x5^2 + 0.5629*x5*x6 - 0.3395*x5*x7 - 0.3277*x5*x8 + 0.06521*x5*x9 + 1.357*x5 - 0.1549*x6^2 - 0.09829*x6*x7 - 0.06071*x6*x8 - 0.2566*x6*x9 - 0.2111*x6 - 0.04051*x7^2 + 0.4806*x7*x8 + 0.2038*x7*x9 + 0.6266*x7 - 0.02868*x8^2 + 0.002664*x8*x9 - 0.1777*x8 - 0.1492*x9^2 - 0.4161*x9 + 30.41;

pLambdaDegree = 0;



% run and verify
lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, [], psy, [], pLambdaDegree, initPhy);
[lpVer, solveResVer, resNorms] = lpVer.solve();

lp4.Lp4Config.printVerifyWithOnlyPsyResult(solveResVer, resNorms);

warning('on')

echo off;
end
