classdef LinearProgram4VerificationBase  < lp4.Lp4AndHlpVerificationBase
    
    properties
        f % 问题1中的f
        
        phy % 问题1中的多项式φ
        lambda % 问题1中的多项式λ
        
        c1Length % 记录约束 1-3 的决策变量 C 的长度
        c2Length
        c3Length
        
        rouIndex = -1
    end % properties
    
    methods
        
        function res = getRouIndex(this)
            res = this.rouIndex;
        end
        
    end % methods
    
    methods (Static)
        
    end % methods (Static)
    
end
