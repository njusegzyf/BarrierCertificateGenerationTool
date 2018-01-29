function [lpVer, solveResVer, resNorms] = runAndVerifyHLPWithIterations2(...
    vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree, pLambdaDegree, pReDegree,...
    maxIterations, initLambdas, initRes)

phys = [];
lambdas = initLambdas;
res = initRes;

for iteration = 1 : maxIterations
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Iteration ', num2str(iteration), ' :']);
    lp4.Lp4Config.displayDelimiterLine();
    
    
    import lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe
    lpVer = HybridLinearProgramVerificationWithGivenLambdaAndRe.createWithRou(...
        vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, degree,...
        lambdas, res);
    
    [lpVer, solveResVer, resNorms] = lpVer.solve();
    
    if solveResVer.hasSolutionWithRou()
        disp('Verify feasible solution succeed, norms :');
        disp(resNorms);
        return;
    elseif ~(solveResVer.hasSolution())
        disp('Unable to find lambda and re for next iterations.');
        return;
    else
        disp(['The rou is: ', num2str(solveResVer.getRou())]);
        lp4.Lp4Config.displayDelimiterLine();
        phys = solveResVer.getPhyExpressions();
    end
    
    
    import lp4.HybridLinearProgramVerificationWithGivenPhy
    lpVer = HybridLinearProgramVerificationWithGivenPhy.createWithRou(...
        vars, stateNum, fs, eps, thetaStateIndex, theta, psys, zetas, guards, pLambdaDegree, pReDegree,...
        phys);
    
    [lpVer, solveResVer, resNorms] = lpVer.solve();
    
    if solveResVer.hasSolutionWithRou()
        disp('Verify feasible solution succeed, norms :');
        disp(resNorms);
        return;
    elseif ~(solveResVer.hasSolution())
        disp('Unable to find lambda and re for next iterations.');
        return;
    else
        disp(['The rou is: ', num2str(solveResVer.getRou())]);
        lp4.Lp4Config.displayDelimiterLine();
        lambdas = solveResVer.getPLmabdaExpressions();
        res = solveResVer.getPReExpressions();
    end
end

lp4.Lp4Config.displayDelimiterLine();
disp(['Over max iteration : ', num2str(maxIterations), ', and no solution is find.']);
lp4.Lp4Config.displayDelimiterLine();

end
