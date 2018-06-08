classdef HybridLinearProgramVerificationWithGivenPhySolveRes < lp4.HybridLinearProgramVerificationSolveResBase
    %HybridLinearProgramVerificationWithGivenPhySolveRes Represents a solve result of a HybridLinearProgramVerificationWithGivenPhy.
    
    methods
        function this = HybridLinearProgramVerificationWithGivenPhySolveRes(linearProgramArg, xArg, fvalArg, exitflagArg, timeArg)
            %HybridLinearProgramVerificationWithGivenPhySolveRes 构造此类的实例
            
            this@lp4.HybridLinearProgramVerificationSolveResBase(linearProgramArg, xArg, fvalArg, exitflagArg, timeArg);
        end
        
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
        
        function printSolution(this)
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
        end
    end
    
    methods (Static = true)
    end
end
