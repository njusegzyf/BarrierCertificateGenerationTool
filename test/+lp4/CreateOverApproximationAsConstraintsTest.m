classdef CreateOverApproximationAsConstraintsTest < matlab.unittest.TestCase

    methods (Test)
	
        function testCreateOverApproximationAsConstraints1(testCase)
            import lp4.createOverApproximationAsConstraints
			import lp4util.Partition
			
			syms p, c, w;
            res = createOverApproximationAsConstraints(p, Partition(0, 4), c, Partition(0, 4), w);
        end

		function testCreateOverApproximation2(testCase)
            import lp4.createOverApproximationAsConstraints
			import lp4util.Partition
			
			syms p, c, w;
            res = createOverApproximationAsConstraints(p, Partition(-1, 1), c, Partition(-1, 1), w);
        end
		
		function testCreateOverApproximation3(testCase)
            import lp4.createOverApproximationAsConstraints
			import lp4util.Partition
			
			syms p, c, w;
            res = createOverApproximationAsConstraints(p, Partition(-4, 4), c, Partition(-4, 4), w);
        end
		
    end % end methods (Test)

    properties
    end
end


