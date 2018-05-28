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
        
        function expr = subAllDecvarValues(this, expr)
            lp = this.linearProgram;
            decvars = lp.decvars;
            decvarsLen = length(decvars);
            decvarValues = this.x;
            
            for i = 1 : decvarsLen
                expr = subs(expr, decvars(i), decvarValues(i));
            end
        end
        
        function expr = subAllIndvarValues(this, expr, indvarValues)
            lp = this.linearProgram;
            indvars = lp.indvars;
            indvarsLen = length(indvars);
            
            for i = 1 : indvarsLen
                expr = subs(expr, indvars(i), indvarValues(i));
            end
        end
        
        function thetaRightExpr = getThetaConstraintRightExpr(this)  
            lp = this.linearProgram;
            thetaRightExpr = lp.exprs(1).polyexpr + lp.phy;
            if lp4.Lp4Config.IS_ADD_EPS_IN_THETA_EXP 
                thetaRightExpr = thetaRightExpr - lp.eps(2);
            end
            
            thetaRightExpr = this.subAllDecvarValues(thetaRightExpr);
        end
        
        function thetaRightExprValue = computeThetaConstraintRightExprValue(this, indvarValues) 
            
            thetaRightExpr = this.getThetaConstraintRightExpr();
            
            thetaRightExpr = this.subAllIndvarValues(thetaRightExpr, indvarValues);
            
            % vpa(x) uses variable-precision floating-point arithmetic (VPA) to evaluate each element of the symbolic input x 
            % to at least d significant digits, where d is the value of the digits function. The default value of digits is 32. 
            thetaRightExprValue = vpa(thetaRightExpr);
        end
        
        function zetaRightExpr = getZetaConstraintRightExpr(this) 
            lp = this.linearProgram;
            lpExprs = lp.exprs;
            % find the zeta constraint by expr name = 'zeta',
            zetaExpr = lpExprs(find(arrayfun(@(x) strcmp(x.name, 'zeta'), lp.exprs), 1));
            zetaRightExpr = zetaExpr.polyexpr - lp.phy - lp.eps(2);
            
            zetaRightExpr = this.subAllDecvarValues(zetaRightExpr);
        end
        
    end
end

