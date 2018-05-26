classdef HybridLinearProgramVerificationBase < lp4.Lp4AndHlpVerificationBase
    
    properties
        stateNum
        thetaStateIndex
        guardNum
        
        phys % 在不同状态下的 phy
        fs % 在不同状态下的 f
        
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
