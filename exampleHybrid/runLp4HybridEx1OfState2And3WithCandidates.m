function [lp, solveRes, lpVer, solveResVer, isVerified, resNorms] = runLp4HybridEx1OfState2And3WithCandidates()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

% get the problem
[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx1OfState2And3Problem();

fixedPhyIndex = 2;
fixedPhy = [1; 1; 1; 1; 1; 1];


% Set the degree of phy, lambda and re
degrees = [2, 3];
pLambdaDegrees = [1]; % [1, 2];
pReDegrees = [1]; % [1, 2];

import lp4util.createRangeCandidates
ranges = [1];
rangeCandicatesCount = length(ranges);
[phyRanges, pLambdaRanges, pReRanges] = createRangeCandidates(ranges, ranges, ranges);

% phys = [1, 0.3, 0.1];
% pLambdas = [1, 0.3, 0.1];
% pRes = [1, 0.3, 0.1];
% rangeCandicatesCount = length(phys);
% [phyRanges, pLambdaRanges, pReRanges] = createRangeCandidates(phys, pLambdas, pRes);



isVerified = false;

for degree = degrees
    for pLambdaDegree = pLambdaDegrees
        for pReDegree = pReDegrees
            
            import lp4.HybridLinearProgram
            lp = HybridLinearProgram.createWithoutRanges(vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards,...
                degree, pLambdaDegree, pReDegree);
            
            for i = 1 : rangeCandicatesCount
                phyRange = phyRanges(i);
                pLambdaRange = pLambdaRanges(i);
                pReRange = pReRanges(i);
                
                lp = lp.setRangeAndWConstraint(phyRange, pLambdaRange, pReRange);
                lp = lp.setPhyExpression(fixedPhyIndex, fixedPhy);
                
                % solve the lp problem
                [lp, solveRes] = lp.solve();
                
                isVerified = false;
                
                % verify the lp problem with the computed lambda
                [lpVer, solveResVer, resNorms] = lp.verifyWithPhy(solveRes);
                import lp4.isResNormsOk
                if (solveRes.hasSolution() && solveResVer.hasSolution() && isResNormsOk(resNorms))
                    isVerified = true;
                    return;
                end
                
                [lpVer, solveResVer, resNorms] = lp.verifyWithLambdaAndRe(solveRes);
                import lp4.isResNormsOk
                if (solveRes.hasSolution() && solveResVer.hasSolution() && isResNormsOk(resNorms))
                    isVerified = true;
                    return;
                end
                
                lp4.Lp4Config.displayDelimiterLine();
                
                % if a large range has no solution, we can skip small ranges
                if ~solveRes.hasSolution()
                    continue;
                end
            end
            
        end
    end
end

warning('on')

echo off;
end
