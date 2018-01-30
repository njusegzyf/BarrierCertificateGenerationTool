classdef CvxSolveRes < lp4util.SolveResBase
    
    properties
    end
    
    methods
        function this = CvxSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime)
            %CvxSolveRes 构造此类的实例
            
            exitflag = lp4util.CvxSolveRes.convertCvxStatusToExitflag(cvxStatus);
            
            %this@lp4util.SolveResBase(linearProgram, x, cvxOptval, exitflag, cvxCpuTime);
            this.linearProgram = linearProgram;
            this.x = x;
            this.fval = cvxOptval;
            this.exitflag = exitflag;
            this.time = cvxCpuTime;
            
            this.output = cvxStatus;
        end
        
    end
    
    methods (Static)
        
        function exitflag = convertCvxStatusToExitflag(cvxStatus)
            % @see http://cvxr.com/cvx/doc/solver.html
            
            if strcmp(cvxStatus, 'Solved') || strcmp(cvxStatus, 'Inaccurate/Solved')
                exitflag = 1;
            elseif strcmp(cvxStatus, 'Unbounded') || strcmp(cvxStatus, 'Infeasible')
                exitflag = -1;
            elseif strcmp(cvxStatus, 'Failed')
                exitflag = -1;
            else
                exitflag = -1;
            end
        end
        
    end
end

