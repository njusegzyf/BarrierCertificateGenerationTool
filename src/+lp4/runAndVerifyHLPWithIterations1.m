function [lpVer, solveResVer, resNorms] = runAndVerifyHLPWithIterations1(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initPhy)

phy = initPhy;
lambda = [];

for iteration = 1 : maxIterations
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Iteration ', num2str(iteration), ' :']);
    lp4.Lp4Config.displayDelimiterLine();
    
    
    import lp4.HybridLinearProgramVerificationWithGivenPhy
    lpVer = HybridLinearProgramVerificationWithGivenPhy.createWithRou(...
        vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree,...
        phy);
    
    [lpVer, solveResVer, resNorms] = lpVer.solve();
    
    if solveResVer.hasSolutionWithRou()
        disp('Verify feasible solution succeed, norms :');
        disp(resNorms);
        return;
    elseif ~(solveResVer.hasSolution())
        disp('Unable to find lambda and re for an next iteration.');
        return;
    else
        disp(['The rou is: ', num2str(solveResVer.getRou())]);
        lp4.Lp4Config.displayDelimiterLine();
        lambda = solveResVer.getPLmabdaExpressions();
        res = solveResVer.getPReExpressions();
    end
    
    
    import lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe
    lpVer = HybridLinearProgramVerificationWithGivenLambdaAndRe.createWithRou(...
        vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree,...
        lambda, res);
    
    [lpVer, solveResVer, resNorms] = lpVer.solve();
    
    if solveResVer.hasSolutionWithRou()
        disp('Verify feasible solution succeed, norms :');
        disp(resNorms);
        return;
    elseif ~(solveResVer.hasSolution())
        disp('Unable to find phy for an next iteration.');
        return;
    else
        disp(['The rou is: ', num2str(solveResVer.getRou())]);
        lp4.Lp4Config.displayDelimiterLine();
        phy = solveResVer.getPhyExpressions();
    end
end

lp4.Lp4Config.displayDelimiterLine();
disp(['Over max iteration : ', num2str(maxIterations), ', and no solution is find.']);
lp4.Lp4Config.displayDelimiterLine();

end
