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
        
    end % methods
    
    methods (Static)
        
    end % methods (Static)
    
end
