classdef LinearProgram4Verification2SolveResult
    %LinearProgram4Verification2SolveResult Represents a solve result of LinearProgram4Verification2.
    
    properties
        linearProgram
        
        % the result of solving the linear programming
        x
        % the value of the objective function fun at the solution x: fval = f'*x.
        fval
        % a value exitflag that describes the exit condition
        exitflag
        % a structure output that contains information about the optimization process
        % output
        
        time
        
        normThreadhold = 0.00001
    end
    
    methods
        function this = LinearProgram4Verification2SolveResult(linearProgramArg, xArg, fvalArg, exitflagArg, timeArg)
            %LinearProgram4Verification2SolveResult 构造此类的实例
            if ~(isa(linearProgramArg, 'LinearProgram4Verification2'))
                error('');
            end
            
            this.linearProgram = linearProgramArg;
            this.x = xArg;
            this.fval = fvalArg;
            this.exitflag = exitflagArg;
            this.time = timeArg;
        end
        
        function res = hasSolution(this)
            res = (this.exitflag == 1);
        end
        
        function res = getPhyCoefficient(this)
            start = this.linearProgram.getPhyCoefficientStart();
            res = this.x(start : start + this.linearProgram.getPhyCoefficientLength() - 1);
        end
        
        function res = getPhyExpression(this)
            lp = this.linearProgram;
            import lp4util.reshapeToVector
            res = reshapeToVector(this.getPhyCoefficient()) * monomials(lp.indvars, 0 : lp.degree);
        end
        
        function res = verify(this)
            for i = 1 : 3
                if ~(this.verifyExpr(i))
                    res = false;
                    return
                end
            end
            
            res = true;
        end
        
        function res = verifyExpr(this, index)
            lp = this.linearProgram;
            expr = lp.exprs(index);
            
            if strcmp(expr.name, 'empty')
                res  = true;
                return;
            end
            
            mid = expr.A * this.x - expr.b;
            res = norm(mid) <= this.normThreadhold;
        end
        
        function res = verifyNorms(this)
            res = [];
            for i = 1 : 3
                res = [res, this.verifyExprNorm(i)];
            end
        end
        
        function res = verifyExprNorm(this, index)
            lp = this.linearProgram;
            expr = lp.exprs(index);
            
            if strcmp(expr.name, 'empty')
                res = 0;
                return;
            end
            
            mid = expr.A * this.x - expr.b;
            res = norm(mid);
        end
        
        function printSolution(this)
            lp = this.linearProgram;
            flag = this.exitflag;
            
            % diaplay code from lp3
            disp('--------------------------------------------------------------');
            disp('The parameter setting:');
            disp(['degree: ', num2str(lp.degree),...
                '; eps1: ',num2str(lp.eps(1)),...
                '; eps2: ',num2str(lp.eps(2))]);
            
            import lp4util.reshapeToVector
            
            if (flag == 1)
                disp('--------------------------------------------------------------');
                disp('The coefficients of function phy is:');
                disp(reshapeToVector(this.getPhyCoefficient()));
                disp('--------------------------------------------------------------');
                disp('The function phy is:');
                disp(this.getPhyExpression());
                disp('--------------------------------------------------------------');
                disp('The computation time is:');
                disp(this.time);
                disp('--------------------------------------------------------------');
            elseif (flag == 0)
                disp('--------------------------------------------------------------');
                disp('Maximum number of iterations reached.');
                disp('--------------------------------------------------------------');
            elseif flag < 0
                disp('--------------------------------------------------------------');
                disp(['The problem with degree ', num2str(lp.degree),' maybe have no solution.']);
                disp('--------------------------------------------------------------');
            else % flag > 1
            end
        end
    end
    
    methods (Static = true)
    end
end

