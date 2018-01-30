classdef HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes
    %HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes Represents a solve result of a HybridLinearProgramVerificationWithGivenLambdaAndRe.
    
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
    end
    
    methods
        function this = HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes(linearProgramArg, xArg, fvalArg, exitflagArg, timeArg)
            %HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes 构造此类的实例
            this.linearProgram = linearProgramArg;
            this.x = xArg;
            this.fval = fvalArg;
            this.exitflag = exitflagArg;
            this.time = timeArg;
        end
        
        function res = hasSolution(this)
            res = (this.exitflag == 1);
        end
        
        function res = hasSolutionWithRou(this)
            res = this.exitflag == 1 && this.getRou() <= 0;
        end
        
        function res = getRou(this)
            if this.linearProgram.isAttachRou
                % as we set linprogF = rou, so we can directly return fval
                res = this.fval;
                % res = this.x(this.linearProgram.decvarsIndexes.rouIndex);
            else
                error('Rou is not used.')
            end
        end
        
        function res = getPhyCoefficient(this, i)
            res = this.x(this.linearProgram.getPhyCoefficientStart(i) : this.linearProgram.getPhyCoefficientEnd(i));
            % res = reshape(res, 1, size(res, 1));
        end
        
        function res = getPhyExpression(this, i)
            lp = this.linearProgram;
            import lp4util.reshapeToVector
            res = reshapeToVector(this.getPhyCoefficient(i)) * monomials(lp.indvars, 0 : lp.degree);
        end
        
        function res = getPhyExpressions(this)
            res(this.linearProgram.stateNum) = this.getPhyExpression(this.linearProgram.stateNum);
            for i = 1 : this.linearProgram.stateNum - 1
                res(i) = this.getPhyExpression(i);
            end
        end
        
        function res = getCGuardCoefficient(this, i)
            decvarIndexes = this.linearProgram.decvarsIndexes;
            res = this.x(decvarIndexes.cGuardStarts(i):  decvarIndexes.cGuardEnds(i));
        end
        
        function guardExpr = getGuardExpr(this, i)
            lp = this.linearProgram;
            guardIndex = 1 + lp.stateNum + i;
            guardExpr = lp.exprs(guardIndex).polyexpr;
            
            for i = 1 : length(lp.decvars)
                decvar = lp.decvars(i);
                decvarValue = this.x(i);
                guardExpr = subs(guardExpr, decvar, decvarValue);
            end
        end
        
        function psyExpr = getPsyExpr(this, i)
            lp = this.linearProgram;
            psyIndex = 1 + i;
            psyExpr = lp.exprs(psyIndex).polyexpr;
            
            for i = 1 : length(lp.decvars)
                decvar = lp.decvars(i);
                decvarValue = this.x(i);
                psyExpr = subs(psyExpr, decvar, decvarValue);
            end
        end
        
        function res = computeExprNorm(this, index)
            lp = this.linearProgram;
            expr = lp.exprs(index);
            
            if strcmp(expr.name, 'empty') || strcmp(expr.type, 'ie')
                res = 0;
                return;
            end
            
            mid = expr.A * this.x - expr.b;
            res = norm(mid);
        end
        
        function res = computeAllExprsNorms(this)
            exprCount = length(this.linearProgram.exprs);
            res(exprCount) = this.computeExprNorm(exprCount);
            for i = 1 : exprCount - 1
                res(i) = this.computeExprNorm(i);
            end
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
                for i = 1 : this.linearProgram.stateNum
                    disp(['The coefficients of function phy', num2str(i), ' is:']);
                    disp(reshapeToVector(this.getPhyCoefficient(i)));
                    disp(['The function phy', num2str(i), ' is:']);
                    disp(reshapeToVector(this.getPhyExpression(i)));
                    disp('');
                end
                
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

