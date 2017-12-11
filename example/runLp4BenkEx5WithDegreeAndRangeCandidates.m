function [lp, solveRes, lpVer, solveResVer, resNorms] = runLp4BenkEx8WithDegreeAndRangeCandidates()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, f, eps, g_theta, g_psy, g_zeta] = getLp4BenkEx8Problem();



% Set the degree of phy and lambda
degrees = [1, 2, 3, 4];
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



% degree: 1; lambda degree: 1
% The function phy is:
% 5542787919849213/9007199254740992 - (4850220904844769*y)/9007199254740992 - (2425110452422385*x)/4503599627370496
% The function lambda is:
% (7883464268065035*x)/9007199254740992 + (7883464268065035*y)/9007199254740992 - 3941244300986739/2251799813685248
 

% --------------------------------------------------------------
% The parameter setting:
% degree: 1; lambda degree: 1; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function phy is:
    % 0.6154   -0.5385   -0.5385

% --------------------------------------------------------------
% The function phy is:
% 5542787919849213/9007199254740992 - (4850220904844769*y)/9007199254740992 - (2425110452422385*x)/4503599627370496
 
% --------------------------------------------------------------
% The coefficients of lambda is:
   % -0.6154         0         0

% --------------------------------------------------------------
% The function lambda is:
% -2771393959924607/4503599627370496
 
% --------------------------------------------------------------
% The rou is:
   % -0.1538

% --------------------------------------------------------------
% The computation time is:
    % 0.3763

% --------------------------------------------------------------
% constraint theta is processed: 2017-12-11 13:02:20
% constraint psy is processed: 2017-12-11 13:02:20
% constraint zeta is processed: 2017-12-11 13:02:20

% No feasible solution found.

% Linprog stopped because no point satisfies the constraints.

% --------------------------------------------------------------
% The parameter setting:
% degree: 1; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The problem with degree 1 maybe have no solution.
% --------------------------------------------------------------
% Verify feasible solution failed.
% constraint theta is processed: 2017-12-11 13:02:28
% constraint psy is processed: 2017-12-11 13:02:29
% constraint zeta is processed: 2017-12-11 13:02:29

% Optimal solution found.

% --------------------------------------------------------------
% The parameter setting:
% lambdaDegree: 1; eps1: 0.0001; eps2: 0.0001
% --------------------------------------------------------------
% The coefficients of function lambda is:
   % -1.7503
    % 0.8752
    % 0.8752

% --------------------------------------------------------------
% The function lambda is:
% (7883464268065035*x)/9007199254740992 + (7883464268065035*y)/9007199254740992 - 3941244300986739/2251799813685248
 
% --------------------------------------------------------------
% The computation time is:
    % 0.1494

% --------------------------------------------------------------
% Verify feasible solution succeed, norms ;
   % 1.0e-15 *

         % 0    0.3513    0.1110
