classdef HybridLinearProgramVerificationSolveResBase < lp4util.SolveResBase
    
    methods
        
        function this = HybridLinearProgramVerificationSolveResBase(linearProgram, x, fval, exitflag, time)
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
            
            if isa(this.linearProgram, 'lp4.HybridLinearProgramVerificationWithGivenLambdaAndRe')
                
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
                
            elseif isa(this.linearProgram, 'lp4.HybridLinearProgramVerificationWithGivenPhy')
                
                lp = this.linearProgram;
                flag = this.exitflag;
                
                % diaplay code from lp3
                disp('--------------------------------------------------------------');
                disp('The parameter setting:');
                disp(['; lambda degree: ', num2str(lp.pLambdaDegree),...
                    '; re degree: ', num2str(lp.pReDegree),...
                    '; eps1: ',num2str(lp.eps(1)),...
                    '; eps2: ',num2str(lp.eps(2))]);
                
                import lp4util.reshapeToVector
                
                if (flag == 1)
                    disp('--------------------------------------------------------------');
                    for i = 1 : this.linearProgram.stateNum
                        disp(['The coefficients of lambda', num2str(i), ' is:']);
                        disp(reshapeToVector(this.getPLmabdaCoefficient(i)));
                        disp(['The function lambda', num2str(i), ' is:']);
                        disp(reshapeToVector(this.getPLmabdaExpression(i)));
                        disp('');
                    end
                    
                    disp('--------------------------------------------------------------');
                    for i = 1 : this.linearProgram.guardNum
                        disp(['The coefficients of re', num2str(i), ' is:']);
                        disp(reshapeToVector(this.getPReCoefficient(i)));
                        disp(['The function re', num2str(i), ' is:']);
                        disp(reshapeToVector(this.getPReExpression(i)));
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
                    disp('The problem maybe have no solution.');
                    disp('--------------------------------------------------------------');
                else % flag > 1
                end
                
            else
                error('Unknown HLP type: %s', class(this.linearProgram));
            end
            
        end
        
    end
    
    methods % for HybridLinearProgramVerificationWithGivenLambdaAndRe
        
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
        
    end % methods HybridLinearProgramVerificationWithGivenLambdaAndRe
    
    methods % for HybridLinearProgramVerificationWithGivenPhy
        
        function res = getPLmabdaCoefficient(this, i)
            res = this.x(this.linearProgram.getPLambdaCoefficientStart(i) : this.linearProgram.getPLambdaCoefficientEnd(i));
            % res = reshape(res, 1, size(res, 1));
        end
        
        function res = getPLmabdaExpression(this, i)
            lp = this.linearProgram;
            import lp4util.reshapeToVector
            res = reshapeToVector(this.getPLmabdaCoefficient(i)) * monomials(lp.indvars, 0 : lp.pLambdaDegree);
        end
        
        function res = getPLmabdaExpressions(this)
            res(this.linearProgram.stateNum) = this.getPLmabdaExpression(this.linearProgram.stateNum);
            for i = 1 : this.linearProgram.stateNum - 1
                res(i) = this.getPLmabdaExpression(i);
            end
        end
        
        function res = getPReCoefficient(this, i)
            res = this.x(this.linearProgram.getPReCoefficientStart(i) : this.linearProgram.getPReCoefficientEnd(i));
            % res = reshape(res, 1, size(res, 1));
        end
        
        function res = getPReExpression(this, i)
            lp = this.linearProgram;
            import lp4util.reshapeToVector
            res = reshapeToVector(this.getPReCoefficient(i)) * monomials(lp.indvars, 0 : lp.pReDegree);
        end
        
        function res = getPReExpressions(this)
            res(this.linearProgram.guardNum) = this.getPReExpression(this.linearProgram.guardNum);
            for i = 1 : this.linearProgram.guardNum - 1
                res(i) = this.getPReExpression(i);
            end
        end
        
    end % methods for HybridLinearProgramVerificationWithGivenPhy
end
