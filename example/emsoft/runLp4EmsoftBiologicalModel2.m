function [lpVer, solveResVer, resNorms] = runLp4EmsoftBiologicalModel2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftBiologicalModel2();

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);
x5 = vars(5);
x6 = vars(6);
x7 = vars(7);
x8 = vars(8);
x9 = vars(9);

degree = 2;
pLambdaDegree = 1;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

% initPhy = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9;
% this is computed from above initPhy with degree = 1
% initPhy = (1446247084809473*x2)/4611686018427387904 + (7231374831680253*x3)/2361183241434822606848 + (940704830558181*x4)/2305843009213693952 + (8317549420363109*x5)/4722366482869645213696 + (1807843820267583*x6)/590295810358705651712 + (8866204539489127*x7)/36893488147419103232 + (2834607139113733*x8)/18446744073709551616 + (8344624299873613*x9)/73786976294838206464 - 2547804529151789/2305843009213693952;

initPhy = -8.967*x1-6.376*x2-10.73*x3-6.947*x4+1.357*x5-.2111*x6+.6266*x7-.1777*x8-.4161*x9+30.41+.5227*x1*x2+2.492*x1*x3+.3089*x1*x4+1.391*x1*x5-.114*x1*x6-.8628*x1*x7+.2594*x1*x8+.1612*x1*x9+.8764*x2*x3+.9352*x2*x4-.3628*x2*x5-0.4419e-1*x2*x6+.4258*x2*x7-0.8721e-1*x2*x8-0.4147e-2*x2*x9+1.022*x3*x4-0.5612e-1*x3*x5-0.2406e-2*x3*x6-1.058*x3*x7+.197*x3*x8+.2396*x3*x9-.1054*x4*x5+.3233*x4*x6+1.207*x4*x7-.4179*x4*x8-.115*x4*x9+.5629*x5*x6-.3395*x5*x7-.3277*x5*x8+0.6521e-1*x5*x9-0.9829e-1*x6*x7-0.6071e-1*x6*x8-.2566*x6*x9+.4806*x7*x8+.2038*x7*x9+0.2664e-2*x8*x9+.4788*x1^2+.5521*x2^2-.1785*x3^2+.2405*x4^2-.4099*x5^2-.1549*x6^2-0.4051e-1*x7^2-0.2868e-1*x8^2-.1492*x9^2;


% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
