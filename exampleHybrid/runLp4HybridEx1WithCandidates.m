function [lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runLp4HybridEx1WithCandidates() 

clear; 
echo on;

% disable warning of `Support of character vectors will be removed in a future release.` 
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx1Problem();


% Set the degree of phy, lambda and re
degrees = [5, 6];
pLambdaDegrees = [1];
pReDegrees = [1];

import lp4util.Partition
phys = [1]; % [1, 0.3, 0.1];
pLambdas = [1]; % [1, 0.3, 0.1];
pRes = [1]; % [1, 0.3, 0.1];
import lp4util.createRangeCandidates
[phyRanges, pLambdaRanges, pReRanges] = createRangeCandidates(phys, pLambdas, pRes);



import lp4.runAndVerifyHLPWithDegreeAndRangeCandidates
[lp, solveRes, lpVer, solveResVer, resNorms, isVerified] = runAndVerifyHLPWithDegreeAndRangeCandidates(...
   vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards,...
    degrees, pLambdaDegrees, pReDegrees,...
    phyRanges, pLambdaRanges, pReRanges);

warning('on')

echo off;
end
