classdef SolveResBase
    
    properties
        linearProgram
        
        % the result of solving the linear programming
        x
        % the value of the objective function fun at the solution x: fval = f'*x.
        fval
        % a value exitflag that describes the exit condition
        exitflag
        % a structure output that contains information about the optimization process output
        output
        
        time
        
        resNorms
        
        normThreadhold = lp4.Lp4Config.RES_NORM_THRESHOLD;
    end
    
    methods
        
        function this = SolveResBase(linearProgram, x, fval, exitflag, time)
            % supports the no argument case with the if statement
            % @see http://cn.mathworks.com/help/matlab/matlab_oop/subclass-constructors.html#brr1h8r
            if nargin ~= 0
                this.linearProgram = linearProgram;
                this.x = x;
                this.fval = fval;
                this.exitflag = exitflag;
                this.time = time;
            end
        end
        
        function res = hasSolution(this)
            res = (this.exitflag == 1);
        end
        
    end
    
    methods (Abstract)
        
        printSolution(this)
        
    end
    
end

