classdef MonomialNumberTest < matlab.unittest.TestCase

    methods (Test)
        function testMonomialNumber1and1(testCase)
            actSolution = monomialNumber(1, 1);
            expSolution = 2;
            testCase.verifyEqual(actSolution, expSolution);
        end
        function testMonomialNumber2and2(testCase)
            actSolution = monomialNumber(2, 2);
            expSolution = 6;
            testCase.verifyEqual(actSolution, expSolution);
        end

    end

    properties
    end
end
