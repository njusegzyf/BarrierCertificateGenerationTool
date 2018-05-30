function [lpVer, solveResVer, resNorms] = verifyLp4EmsoftInvertedPendulumLinearizedWithOnlyPsy()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, f, eps, ~, psy, ~] = getLp4EmsoftInvertedPendulumLinearizedProblem();

x1 = vars(1);
x2 = vars(2);
x3 = vars(3);
x4 = vars(4);

initPhy = 2.6132091261160179662681457557483*x1^2 + 5.3337517151510649426882082480006*x1*x2 - 26.429614724711921525113211828284*x1*x3 - 6.6206801732729632092855354130734*x1*x4 + 4.5817259017870872739308651944157*x2^2 - 48.233137865255159226762771140784*x2*x3 - 12.23395703407775769733234483283*x2*x4 + 139.29609345636919215394300408661*x3^2 + 71.639069391630783911750768311322*x3*x4 + 9.4371334500327161975974377128296*x4^2 - 0.1;

pLambdaDegree = 0;



% run and verify
lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, [], psy, [], pLambdaDegree, initPhy);
[lpVer, solveResVer, resNorms] = lpVer.solve();

lp4.Lp4Config.printVerifyWithOnlyPsyResult(solveResVer, resNorms);

warning('on')

echo off;
end
