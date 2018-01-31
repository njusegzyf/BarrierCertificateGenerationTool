function [lpVer, solveResVer, resNorms] = runAndVerifyLp4WithIterations1(...
    vars, f, eps, theta, psy, zeta, degree, pLambdaDegree,...
    maxIterations, initPhy)

phy = initPhy;
lambda = 0;

for iteration = 1 : maxIterations
    lp4.Lp4Config.displayDelimiterLine();
    disp(['Iteration ', num2str(iteration), ' :']);
    lp4.Lp4Config.displayDelimiterLine();
    
    
    lpVer = lp4.LinearProgram4Verification3.createWithRou(vars, f, eps, theta, psy, zeta, pLambdaDegree, phy);
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
        lambda = solveResVer.getLambdaExpression();
    end
    
    
    lpVer = lp4.LinearProgram4Verification2.createWithRou(vars, f, eps, theta, psy, zeta, degree, lambda, 0);
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
        phy = solveResVer.getPhyExpression();
    end
end

lp4.Lp4Config.displayDelimiterLine();
disp(['Over max iteration : ', num2str(maxIterations), ', and no solution is find.']);
lp4.Lp4Config.displayDelimiterLine();

end
