function [lpVer, solveResVer, resNorms] = runLp4HybridEx4CvxVerificationWithGivenPhy()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx4Problem();
x1 = vars(1);
x2 = vars(2);
x3 = vars(3);

% Set the degree
pLambdaDegree = 1;
pReDegree = 1;

phys = [(7486591293511275*x1)/4611686018427387904 - (2445627490774475*x2)/36893488147419103232 - (5240093099737441*x3)/590295810358705651712 + (932379137192507*x1*x2)/18446744073709551616 + (8093276349972299*x1*x3)/4722366482869645213696 - (2829712385172479*x2*x3)/9444732965739290427392 - (2700021363492765*x1^2)/18446744073709551616 + (7185206497060383*x2^2)/295147905179352825856 + (3578897891431351*x3^2)/75557863725914323419136 + 7817122460283389/1152921504606846976,...
    (1469940654599107*x1*x2)/590295810358705651712 - (5240093099737441*x3)/590295810358705651712 - (5240093099737441*x2)/295147905179352825856 + (4046638174986145*x1*x3)/2361183241434822606848 - (5659424770344957*x2*x3)/18889465931478580854784 - (1741895000713781*x1^2)/1180591620717411303424 - (4236832819068933*x2^2)/2361183241434822606848 - (6024342541130967*x3^2)/18889465931478580854784 - 7978039774514801/9223372036854775808];

% Note: If we use matlab built-in `linprog` in Matlab 2017b, we will get an error of:
% ¥ÌŒÛ π”√ linprog (line 336)
% LINPROG encountered an internal error that caused the solution to lose feasibility. We are sorry for the
% inconvenience.
% Please contact technical support for assistance with your problem, quoting the code "-2_4".

lpVer = lp4.HybridLinearProgramVerificationWithGivenPhy.createWithRou(vars, stateNum,...
    fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree,...
    phys);

[lpVer, solveResVer, resNorms] = lpVer.solveWithCvx();
lp4.Lp4Config.displaySolveRes(solveResVer, resNorms);

end
