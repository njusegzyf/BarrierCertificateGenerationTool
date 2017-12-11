classdef MonomialNumberTest < matlab.unittest.TestCase
    
    methods (Test)
        
        function testMonomialNumber1and1(testCase)
            actualNumber = monomialNumber(1, 1);
            testCase.verifyEqual(actualNumber, 2);
        end
        
        function testMonomialNumber3and1(testCase)
            actualNumber = monomialNumber(3, 1);
            testCase.verifyEqual(actualNumber, 4);
        end
        
        function testMonomialNumber2and2(testCase)
            actualNumber = monomialNumber(2, 2);
            testCase.verifyEqual(actualNumber, 6);
        end
        
    end % methods (Test)
    
end
