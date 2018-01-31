function [lpVer, solveResVer, resNorms] = runAndVerifyLp4WithRounds(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...% the problem
    lambdaStart, lambdaEnd,... % the range of initLambda
    randomStartCount, randomStartBeginIndex, maxIterations, randomSeed,...
    isStopOnError)

% generate random ints that are random seeds for generate random init expression
rng(randomSeed);
randomSeeds = randi(10000000, 1, randomStartCount, 'int32');

for i = randomStartBeginIndex : randomStartCount
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Begin with random start ', num2str(i), ' :']);
    lp4.Lp4Config.displayDelimiterLine();
    
    import lp4util.generateRamdomExprs
    randomExprs = generateRamdomExprs(vars, pLambdaDegree, 1, lambdaStart, lambdaEnd, randomSeeds(1, i));
    initLambda = randomExprs(1); 
    
    
    try
        [lpVer, solveResVer, resNorms] = runAndVerifyLp4WithIterations2(...
            vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
            maxIterations, initLambda);
        
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
