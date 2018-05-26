classdef HybridLinearProgramVerificationBase < lp4.Lp4AndHlpVerificationBase
    
    properties
        stateNum
        thetaStateIndex
        guardNum
        
        phys % �ڲ�ͬ״̬�µ� phy
        fs % �ڲ�ͬ״̬�µ� f
        
        psys
        zetas
        guards
        
        decvarsIndexes
    end % properties
    
    methods
        
        % override
        function  res = getRouIndex(this)
            res = this.decvarsIndexes.rouIndex;
        end
        
        % override
        function res = getCStart(this)
            res = this.decvarsIndexes.cThetaStart(1);
        end
        
        function res = nextExprNumIndex(this)
            res = length(this.exprs) + 1;
        end
        
        function cvxSolveRes = createCvxSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime)
            cvxSolveRes = lp4.HybridLinearProgramCvxVerificationSolveRes(linearProgram, x, cvxOptval, cvxStatus, cvxCpuTime);
        end
        
    end % methods
    
    methods (Static)
        
    end % methods (Static)
    
end
