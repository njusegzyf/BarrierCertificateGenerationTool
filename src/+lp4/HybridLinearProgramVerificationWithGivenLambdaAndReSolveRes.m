classdef HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes < lp4.HybridLinearProgramVerificationSolveResBase
    %HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes Represents a solve result of a HybridLinearProgramVerificationWithGivenLambdaAndRe.
    
    methods
        
        function this = HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes(linearProgramArg, xArg, fvalArg, exitflagArg, timeArg)
            %HybridLinearProgramVerificationWithGivenLambdaAndReSolveRes 构造此类的实例
            
            this@lp4.HybridLinearProgramVerificationSolveResBase(linearProgramArg, xArg, fvalArg, exitflagArg, timeArg);
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

