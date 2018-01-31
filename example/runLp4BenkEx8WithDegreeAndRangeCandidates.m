function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4BenkEx8WithDegreeAndRangeCandidates()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx8Problem();



% Set the degree of phy and lambda
degrees = [2, 3, 4];
pLambdaDegrees = [1, 2, 3];

ranges = [1, 0.5, 0.3, 0.15, 0.1];
import lp4util.createRangeCandidates
[phyRanges, pLambdaRanges, phyRangesInVerify] = createRangeCandidates(ranges, ranges, 0);

% run and verify
import lp4.runAndVerifyWithDegreeAndRangeCandidates
[lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyWithDegreeAndRangeCandidates(...
    vars, f, eps, g_theta, g_psy, g_zeta,...
    degrees, pLambdaDegrees,...
    phyRanges, pLambdaRanges, phyRangesInVerify);

warning('on')

echo off;
end



% use phy = [(x + 1)/2, (y + 1)/2], zeta = [2*x - 1, 2*y - 1]
%
% The parameter setting:
% degree: 2; lambda degree: 1; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function phy is:
%     0.6739   -0.7910   -0.7304   -0.4579   -0.7569   -0.4579
% 
% --------------------------------------------------------------
% The function phy is:
% 3034755492373117/4503599627370496 - (6579169656844065*y)/9007199254740992 - (3408835521268723*x*y)/4503599627370496 - (8249129716660465*x^2)/18014398509481984 - (515570607291279*y^2)/1125899906842624 - (7125087829866755*x)/9007199254740992
%  
% --------------------------------------------------------------
% The coefficients of lambda is:
%    -0.6739    0.4579    0.4579
% 
% --------------------------------------------------------------
% The function lambda is:
% (4124564858330233*x)/9007199254740992 + (8249129716660463*y)/18014398509481984 - 3034755492373117/4503599627370496
%  
% --------------------------------------------------------------
% The rou is:
%    -0.0842
% 
% --------------------------------------------------------------
% The computation time is:
%     0.0511
% 
% --------------------------------------------------------------
% constraint theta is processed: 2017-12-11 15:22:08
% constraint psy is processed: 2017-12-11 15:22:09
% constraint zeta is processed: 2017-12-11 15:22:10
% 
% Optimal solution found.
% 
% --------------------------------------------------------------
% The parameter setting:
% degree: 2; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function phy is:
%    1.0e-03 *
% 
%     0.5463   -0.6706   -0.6675   -0.0293    0.1459   -0.0257
% 
% --------------------------------------------------------------
% The function phy is:
% (1346151138084367*x*y)/9223372036854775808 - (1539103407357769*y)/2305843009213693952 - (6184794268706823*x)/9223372036854775808 - (8634904098623183*x^2)/295147905179352825856 - (7578622606194673*y^2)/295147905179352825856 + 5038397138243771/9223372036854775808
%  
% --------------------------------------------------------------
% The computation time is:
%     0.0737
% 
% --------------------------------------------------------------
% Verify feasible solution succeed, norms ;
%    1.0e-17 *
% 
%     0.3258    0.0455    0.0122
