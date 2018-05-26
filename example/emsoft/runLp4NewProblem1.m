function [lpVer, solveResVer, resNorms] = runLp4NewProblem1()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftNewProblem1();

x = vars(1);
y = vars(2);
z = vars(3);
v = vars(4);

degree = 2;
pLambdaDegree = 1;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

% initPhy = x + y + z + v;

% the following phy is computed by the above random init phy
initPhy = (2744160047761963*v)/35184372088832 + (328848242599795*x)/4398046511104 + (1285563093942491*y)/35184372088832 + (7348710066222461*z)/70368744177664 - (3902431151190319*v*x)/281474976710656 + (4397682911166647*v*y)/9007199254740992 - (764609534775031*v*z)/70368744177664 - (8035690325914565*x*y)/2251799813685248 - (1811114384713785*x*z)/70368744177664 - (4818087011388709*y*z)/4503599627370496 - (6074144109018725*v^2)/9007199254740992 - (3805507662146317*x^2)/562949953421312 + (8686983726165971*y^2)/18014398509481984 - (1962336408668881*z^2)/562949953421312 + 5414103506138829/17592186044416;
initPhy = initPhy / 1e2;


% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
