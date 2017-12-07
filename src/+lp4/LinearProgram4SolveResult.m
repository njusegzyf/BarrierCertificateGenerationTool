classdef LinearProgram4SolveResult
    %LinearProgram4SolveResult Represents a solve result of a LinearProgram4.
    
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
    end
    
    methods
        function this = LinearProgram4SolveResult(linearProgramArg, xArg, fvalArg, exitflagArg, timeArg)
            %LINEARPROGRAM4SOLUTION 构造此类的实例
            this.linearProgram = linearProgramArg;
            this.x = xArg;
            this.fval = fvalArg;
            this.exitflag = exitflagArg;
            this.time = timeArg;
        end
        
        function res = hasSolution(this)
            res = (this.exitflag == 1) && (this.getRou <= 0);
        end
        
        function res = getPhyCoefficient(this)
            start = this.linearProgram.getPhyCoefficientStart();
            res = this.x(start : start + this.linearProgram.getPhyCoefficientLength() - 1);
            % res = reshape(res, 1, size(res, 1));
        end
        
        function res = getPhyExpression(this)
            lp = this.linearProgram;
            import lp4util.reshapeToVector
            res = reshapeToVector(this.getPhyCoefficient()) * monomials(lp.indvars, 0 : lp.degree);
        end
        
        function res = getRou(this)
            res = this.x(this.linearProgram.getRouIndex());
        end
        
        function res = getPLmabdaCoefficient(this)
            start = this.linearProgram.getPLambdaCoefficientStart();
            res = this.x(start : start + this.linearProgram.getPLambdaCoefficientLength() - 1);
            % res = reshape(res, 1, size(res, 1));
        end
        
        function res = getPLmabdaExpression(this)
            lp = this.linearProgram;
            import lp4util.reshapeToVector
            res = reshapeToVector(this.getPLmabdaCoefficient()) * monomials(lp.indvars, 0 : lp.pLambdaDegree);
        end
        
        function res = getW(this)
            start = this.linearProgram.getWStart();
            res = this.x(start : start + this.linearProgram.getWLength() - 1);
            % res = reshape(res, 1, size(res, 1));
        end
        
        function printSolution(this)
            lp = this.linearProgram;
            flag = this.exitflag;
            
            % diaplay code from lp3
            disp('--------------------------------------------------------------');
            disp('The parameter setting:');
            disp(['degree: ', num2str(lp.degree),...
                '; lambda degree: ', num2str(lp.pLambdaDegree),...
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
                disp('The coefficients of lambda is:');
                disp(reshapeToVector(this.getPLmabdaCoefficient()));
                disp('--------------------------------------------------------------');
                disp('The function lambda is:');
                disp(this.getPLmabdaExpression());
                disp('--------------------------------------------------------------');
                %                 disp('The w is:');
                %                 disp(this.getW());
                %                 disp('--------------------------------------------------------------');
                disp('The rou is:');
                disp(this.getRou());
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

