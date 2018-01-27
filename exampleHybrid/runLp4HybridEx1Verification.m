function [lpVer, solveResVer] = runLp4HybridEx1Verification()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards] = getLp4HybridEx1Problem();
t = vars(1);
T = vars(2);

% Set the degree of lambda and re
pLambdaDegree = 4;
pReDegree = 4;

phys = [- (7464992003993547*T)/1152921504606846976 - (3*t)/4 - (1190028783492395*T*t)/18014398509481984 - (8308165003744469*T^2)/72057594037927936 - (2375838052912221*t^2)/36028797018963968 - 3/4,...
    (2702247835339475*T^2)/18014398509481984 + (3*T*t)/4 - (6651935653501687*T)/72057594037927936 + (3737*t^2)/6000 - (1487*t)/3000 + 7/40,...
    (4572219731706301*T^2)/72057594037927936 - (5143747198169589*T)/18014398509481984 - (6755421959053885*t)/9007199254740992 - (2951479051793529*t^2)/1180591620717411303424 - 1501/10000];



import lp4.HybridLinearProgramVerificationWithGivenPhy
lpVer = HybridLinearProgramVerificationWithGivenPhy.create(vars, stateNum,...
    fs, eps, thetaStateIndex, zetaStateIndex, theta, psys, zeta, guards, pLambdaDegree, pReDegree,...
    phys);

[lpVer, solveResVer] = lpVer.solve();



warning('on')
echo off;

end
