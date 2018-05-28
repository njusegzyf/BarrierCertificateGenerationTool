function [lpVer, solveResVer, resNorms] = runLp4EmsoftNewProblem2()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftNewProblem2();

x = vars(1);
y = vars(2);
z = vars(3);
v = vars(4);

degree = 2;
pLambdaDegree = 1;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

% initPhy = x + y + z + v;

% the following phy is computed by the above random init phy
initPhy = (1176760725788843*v)/2199023255552 + (3272618591852713*x)/2199023255552 + (3234017500940035*y)/281474976710656 - (643431063311573*z)/1099511627776 - (1016693353688515*v*x)/17592186044416 - (1534047749038879*v*y)/562949953421312 - (7347958985083897*v*z)/70368744177664 + (7486733131878591*x*y)/2251799813685248 - (2282671951936441*x*z)/4398046511104 + (6583499454999401*y*z)/288230376151711744 - (5765905532721907*v^2)/140737488355328 - (966880209303797*x^2)/2199023255552 + (7971840189206019*y^2)/4503599627370496 + (4193319763494715*z^2)/70368744177664 + 5094347434623321/1099511627776;
initPhy = initPhy / 1e3;



% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
