function [lpVer, solveResVer, resNorms] = runLp4EmsoftNewProblem3()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, theta, psy, zeta] = getLp4EmsoftNewProblem3();

x = vars(1);
y = vars(2);
z = vars(3);
v = vars(4);
w = vars(5);

degree = 2;
pLambdaDegree = 1;

maxIterations = lp4.Lp4Config.DEFAULT_MAX_ITERATION_COUNT;

% initPhy = x + y + z + v + w;

% the following phy is computed by the above random init phy
initPhy = (4760192614766915*v)/17592186044416 + (8562038396685589*w)/17592186044416 + (5074643355880943*x)/35184372088832 + (7349978030594931*y)/17592186044416 + (6338322358467341*z)/1125899906842624 - (3902071258851879*v*w)/140737488355328 - (5225173670526585*v*x)/281474976710656 - (4444938175733429*v*y)/140737488355328 + (8047516919759979*w*x)/562949953421312 - (6132146548158533*v*z)/2251799813685248 - (5659552213311175*w*y)/140737488355328 - (5776084054539495*w*z)/70368744177664 - (6138819489580131*x*y)/281474976710656 + (3556350183066907*x*z)/562949953421312 - (659959958999623*y*z)/281474976710656 + (4731584590749741*v^2)/562949953421312 - (125374736361343*w^2)/17592186044416 - (1262207785233091*x^2)/140737488355328 - (1508045907475231*y^2)/70368744177664 + (1815203345525667*z^2)/4503599627370496 + 337451317771959/137438953472;
initPhy = initPhy / 1e3;


% run and verify
[lpVer, solveResVer, resNorms] = lp4.runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy);

warning('on')

echo off;
end
