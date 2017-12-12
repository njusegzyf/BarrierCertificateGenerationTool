classdef LinearProgram4Verification2Test < matlab.unittest.TestCase
    
    methods (Test)
        
        function testExample7(testCase)
            
            [vars, f, eps, g_theta, g_psy, g_zeta] = getLp4Example7Problem();
            
            phyDegree = 2;
            lambda = 0;

            import lp4.LinearProgram4Verification2
            lp = LinearProgram4Verification2(vars);
            
            lp.f = f;
            lp.eps = eps;
            
            % Set the degree of phy
            lp = lp.setDegreeAndInit(phyDegree);
            
            lp.lambda = lambda;
            
            lp = lp.setThetaConstraint(g_theta);
            lp = lp.setPsyConstraint(g_psy);
            lp = lp.setZetaConstraint(g_zeta);
            lp = lp.generateEqsForConstraint1To3();
            
            lp = lp.setDevVarsConstraint();
            
            % solve the lp problem
            [lp, solveRes, resNorms] = lp.solve();
            
            testCase.verifyTrue(solveRes.hasSolution());
        end
        
    end % methods (Test)
end

