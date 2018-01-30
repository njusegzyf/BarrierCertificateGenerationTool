function [lpVer, solveResVer, resNorms] = runLp4HybridEx1WithRounds()

clear;
echo on;

% disable warning of `Support of character vectors will be removed in a future release.`
% which is produced by function `monomials`.
warning('off')

[vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards] = getLp4HybridEx1Problem();

% Set the degree
degree = 4;
pLambdaDegree = 2;
pReDegree = 2;

randomStartCount = 100;
randomStartBeginIndex = 1;
maxIterations = 20;

randomSeed = 0;
rng(randomSeed);
randomSeeds = rand(2, randomStartCount);

lambdaStart = -1;
lambdaEnd = 1;
resStart = 0;
resEnd = 1;

for i = randomStartBeginIndex : randomStartCount
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Begin with random start ', num2str(i), ' :']);
    lp4.Lp4Config.displayDelimiterLine();
    
    import lp4util.generateRamdomExprs
    initLambdas = generateRamdomExprs(vars, pLambdaDegree, stateNum, lambdaStart, lambdaEnd, randomSeeds(1, i));
    initRes = generateRamdomExprs(vars, degree, length(guards), resStart, resEnd, randomSeeds(2, i));
    
    try
        [lpVer, solveResVer, resNorms] = lp4.runAndVerifyHLPWithIterations2(...
            vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
            maxIterations, initLambdas, initRes);
        
        if solveResVer.hasSolutionWithRou()
            return;
        end
    catch
        % continue if an error occurs in one iteration
    end
end

warning('on')
echo off;

end
