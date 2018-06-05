classdef Partition
    %PARTITION Represents a partition.
    
    properties
        boundLow
        boundHigh
    end
    
    methods
        % Constructor
        function this = Partition(boundLowArg, boundHighArg)
            this.boundLow = boundLowArg;
            this.boundHigh = boundHighArg;
        end
        
        function res = splitPartition(this)
            boundMiddle = (this.boundLow + this.boundHigh) / 2;
            import lp4util.Partition
            res = [Partition(this.boundLow, boundMiddle), Partition(boundMiddle, this.boundHigh)];
        end
        
        function constraints = createConstraints(this, p)
            constraint1 = this.boundLow - p;
            constraint2 = p - this.boundHigh;
            constraints = [constraint1, constraint2];
        end
        
        function expr = createConstraintsAndAddToExpr(this, p, decvars, expr)
            constraints = this.createConstraints(p);
            % constraints = arrayfun(@(x) x == 0, constraints);
            [ A, b ] = equationsToMatrix(constraints, decvars);
            % Note: `A` and `b` is of type `sym`, so we convert them to `double`
            A = double(A);
            b = double(b);
            expr.A = [expr.A; A];
            expr.b = [expr.b; b];
        end
        
        function expr = createConstraintsWithRouAndAddToExpr(this, p, decvars, rouIndex, expr)
            constraints = this.createConstraints(p);
            % constraints = arrayfun(@(x) x == 0, constraints);
            [ A, b ] = equationsToMatrix(constraints, decvars);
            % Note: `A` and `b` is of type `sym`, so we convert them to `double`
            A = double(A);
            b = double(b);
            
            aLength = size(A, 1);
            % add rou to the left part of each constraint
            for i = 1 : 1 : aLength
               A(i, rouIndex) = -1;
            end
            
            expr.A = [expr.A; A];
            expr.b = [expr.b; b];
        end
        
        function res = middlePoint(this)
            res = (this.boundLow + this.boundHigh)/2;
        end
        
    end % end methods
    
    methods (Static = true)
        
        function res = createPartitionFromMidValue(mid, delta)
            if delta < 0
                error('');
            end
            
            import lp4util.Partition
            res = Partition(mid - delta, mid + delta);
        end
        
    end % end methods (Static = true) 
    
end
