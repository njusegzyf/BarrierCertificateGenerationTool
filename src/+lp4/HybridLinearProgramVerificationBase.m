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
        
    end % methods
    
    methods (Static)
        
    end % methods (Static)
    
end
