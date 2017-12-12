classdef LinearProgram4Verification3SolveResult
    %LinearProgram4Verification3SolveResult Represents a solve result of LinearProgram4Verification3.
    
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
        function this = LinearProgram4Verification3SolveResult(linearProgramArg, xArg, fvalArg, exitflagArg, timeArg)
            %LinearProgram4Verification2SolveResult 构造此类的实例
            if ~(isa(linearProgramArg, 'LinearProgram4Verification3'))
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
        
        function res = getLambdaCoefficient(this)
            start = this.linearProgram.getLambdaCoefficientStart();
            res = this.x(start : start + this.linearProgram.getLambdaCoefficientLength() - 1);
        end
        
        function res = getLambdaExpression(this)
            lp = this.linearProgram;
            import lp4util.reshapeToVector
            res = reshapeToVector(this.getLambdaCoefficient()) * monomials(lp.indvars, 0 : lp.lambdaDegree);
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
            
            mid = expr.A * this.x - expr.b;
            res = norm(mid);
        end
        
        function printSolution(this)
            lp = this.linearProgram;
            flag = this.exitflag;
            
            if this.hasSolution()
                disp('Verification succeeded.');
            else
                disp('Verification failed.');
            end
            disp(['The parameter setting:',...
                'lambda degree: ', num2str(lp.lambdaDegree),...
                '; eps1: ',num2str(lp.eps(1)),...
                '; eps2: ',num2str(lp.eps(2))]);
            
            import lp4.Lp4Config
            if this.hasSolution() || Lp4Config.IS_PRINT_FAILED_VERIFICATION_INFO
                % diaplay code from lp3
                disp('--------------------------------------------------------------');
                disp('The parameter setting:');
                disp(['lambdaDegree: ', num2str(lp.lambdaDegree),...
                    '; eps1: ',num2str(lp.eps(1)),...
                    '; eps2: ',num2str(lp.eps(2))]);
                
                if (flag == 1)
                    disp('--------------------------------------------------------------');
                    disp('The coefficients of function lambda is:');
                    disp(this.getLambdaCoefficient());
                    disp('--------------------------------------------------------------');
                    disp('The function lambda is:');
                    disp(this.getLambdaExpression());
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
                    disp(['The problem with lambdaDegree ', num2str(lp.lambdaDegree),' maybe have no solution.']);
                    disp('--------------------------------------------------------------');
                else % flag > 1
                end
            end
        end % function printSolution
        
    end % methods
    
    methods (Static = true)
    end
end

