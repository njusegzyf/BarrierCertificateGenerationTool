classdef Lp4AndHlpVerificationSolveResBase < lp4util.SolveResBase
    
    methods
        
        function this = Lp4AndHlpVerificationSolveResBase(linearProgram, x, fval, exitflag, time)
            % supports the no argument case with the if statement
            % @see http://cn.mathworks.com/help/matlab/matlab_oop/subclass-constructors.html#brr1h8r

            % since we can not call an super class constructor in if statement, 
            % so we always call no argument super class constructor and then init `this` if we are not in no argument case.  
            this@lp4util.SolveResBase();
            if nargin ~= 0
                this.linearProgram = linearProgram;
                this.x = x;
                this.fval = fval;
                this.exitflag = exitflag;
                this.time = time;
            end
        end
        
        function res = hasSolutionWithRou(this)
            res = this.exitflag == 1 && this.getRou() <= lp4.Lp4Config.ROU_THRESHOLD; % this.getRou() <= 0;
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
        
        function expr = computeEq1Right(this, indvarValues) 
            lp = this.linearProgram;
            expr = lp.exprs(1).polyexpr + lp.phy;
            
            decvarLen = length(lp.decvars);
            for i = 1 : decvarLen
                decvar = lp.decvars(i);
                decvarValue = this.x(i);
                expr = subs(expr, decvar, decvarValue);
            end
            
            indvarLen = length(lp.indvars);
            for i = 1 : indvarLen
                indvar = lp.indvars(i);
                indvarValue = indvarValues(i);
                expr = subs(expr, indvar, indvarValue);
            end
            
            expr = vpa(expr);
        end
        
    end
end

