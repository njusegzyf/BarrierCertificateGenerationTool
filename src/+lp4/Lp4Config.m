classdef Lp4Config
    %LP4CONFIG Stores the configurations.
    
    properties (Constant)
        IS_DEBUG = true;
        IS_VERBOSE = false;

        VERIFICATION_C_DEGREE_INC = 2;
    end
    
    methods (Static)
        function res = isDebug()
            import lp4.Lp4Config
            res = Lp4Config.IS_DEBUG;
        end
        
        function res = isVerbose()
            import lp4.Lp4Config
            res = Lp4Config.IS_VERBOSE;
        end
    end
end

