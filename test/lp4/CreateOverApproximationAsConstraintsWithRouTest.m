classdef CreateOverApproximationAsConstraintsWithRouTest < matlab.unittest.TestCase
    
    methods (Test)
        
        function testCreation1(testCase)
            import lp4.createOverApproximationAsConstraintsWithRou
            import lp4util.Partition
            
            syms p c w rou;
            res = createOverApproximationAsConstraintsWithRou(p, Partition(0, 4), c, Partition(0, 4), w, rou);
        end
        
        function testCreation2(testCase)
            import lp4.createOverApproximationAsConstraintsWithRou
            import lp4util.Partition
            
            syms p c w rou;
            res = createOverApproximationAsConstraintsWithRou(p, Partition(-1, 1), c, Partition(-1, 1), w, rou);
        end
        
        function testCreation3(testCase)
            import lp4.createOverApproximationAsConstraintsWithRou
            import lp4util.Partition
            
            syms p c w rou;
            res = createOverApproximationAsConstraintsWithRou(p, Partition(-4, 4), c, Partition(-4, 4), w, rou);
        end
        
        function testValueInConstraints1(testCase)
            import lp4.createOverApproximationAsConstraintsWithRou
            import lp4util.Partition
            
            syms p c w rou;
            decvars = [p, c, w, rou];
            constraints = createOverApproximationAsConstraintsWithRou(p, Partition(-4, 4), c, Partition(-4, 4), w, rou);
            [constraintsIeqA, constraintsIeqb] = equationsToMatrix(constraints, decvars);
            constraintsIeqA = double(constraintsIeqA);
            constraintsIeqb = double(constraintsIeqb);
            
            eqA = ...
                [1, 0, 0, 0;
                0, 1, 0, 0;
                0, 0, 1, 0];
            
            f = zeros(1, 4);
            
            % test p = 0, c = 0, w = 0
            eqb = [0; 0; 0];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
            
            % test p = -4, c = -4, w = 16
            eqb = [-4; -4; 16];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
            
            % test p = 4, c = 4, w = 16
            eqb = [4; 4; 16];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
            
            % test p = 4, c = -4, w = -16
            eqb = [4; -4; -16];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
            
            % test p = 2, c = 2, w = 4
            eqb = [2; 2; 4];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
            
            % test p = -2, c = 2, w = -4
            eqb = [-2; 2; -4];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
        end
        
        function testValueInConstraints2(testCase)
            import lp4.createOverApproximationAsConstraintsWithRou
            import lp4util.Partition
            
            syms p c w rou;
            decvars = [p, c, w, rou];
            constraints = createOverApproximationAsConstraintsWithRou(p, Partition(0, 4), c, Partition(-4, 0), w, rou);
            [constraintsIeqA, constraintsIeqb] = equationsToMatrix(constraints, decvars);
            constraintsIeqA = double(constraintsIeqA);
            constraintsIeqb = double(constraintsIeqb);
            
            eqA = ...
                [1, 0, 0, 0;
                0, 1, 0, 0;
                0, 0, 1, 0];
            
            f = zeros(1, 4);
            
            % test p = 0, c = 0, w = 0
            eqb = [0; 0; 0];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
            
            % test p = 4, c = -4, w = -16
            eqb = [4; -4; -16];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
 
            % test p = 2, c = -2, w = -4
            eqb = [2; -2; -4];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
        end
        
         function testValueInConstraints3(testCase)
            import lp4.createOverApproximationAsConstraintsWithRou
            import lp4util.Partition
            
            syms p c w rou;
            decvars = [p, c, w, rou];
            constraints = createOverApproximationAsConstraintsWithRou(p, Partition(-4, 0), c, Partition(-4, 0), w, rou);
            [constraintsIeqA, constraintsIeqb] = equationsToMatrix(constraints, decvars);
            constraintsIeqA = double(constraintsIeqA);
            constraintsIeqb = double(constraintsIeqb);
            
            eqA = ...
                [1, 0, 0, 0;
                0, 1, 0, 0;
                0, 0, 1, 0];
            
            f = zeros(1, 4);
            
            % test p = 0, c = 0, w = 0
            eqb = [0; 0; 0];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
            
            % test p = -4, c = -4, w = 16
            eqb = [-4; -4; 16];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
 
            % test p = -2, c = -2, w = 4
            eqb = [-2; -2; 4];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
            
            % test p = -4, c = 0, w = 0
            eqb = [-4; 0; 0];
            [~, ~, exitflag] = linprog(f, constraintsIeqA, constraintsIeqb, eqA, eqb);
            testCase.verifyEqual(exitflag, 1);
        end
        
    end % end methods (Test)
    
end
