function [lpVer, solveResVer, resNorms] = runAndVerifyHLPWithRounds(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards,... % the problem
    degree, pLambdaDegree, pReDegree,...
    lambdaStart, lambdaEnd, resStart, resEnd,... % the range of initLambda and initRe
    randomStartCount, randomStartBeginIndex, maxIterations, randomSeed,...
    isStopOnError)

% generate random ints that are random seeds for generate random init expression
rng(randomSeed);
randomSeeds = randi(10000000, 2, randomStartCount, 'int32');

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
    catch err
        disp(err);
        if isStopOnError
            return
        else
            % continue if an error occurs in one iteration
        end
    end
end

end
