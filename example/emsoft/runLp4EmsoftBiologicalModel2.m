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
initPhy = (1446247084809473*x2)/4611686018427387904 + (7231374831680253*x3)/2361183241434822606848 + (940704830558181*x4)/2305843009213693952 + (8317549420363109*x5)/4722366482869645213696 + (1807843820267583*x6)/590295810358705651712 + (8866204539489127*x7)/36893488147419103232 + (2834607139113733*x8)/18446744073709551616 + (8344624299873613*x9)/73786976294838206464 - 2547804529151789/2305843009213693952;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
