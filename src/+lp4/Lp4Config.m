classdef Lp4Config
    %LP4CONFIG Stores the configurations.
    
    properties (Constant)
        IS_DEBUG = true;
        IS_VERBOSE = false;
		
		RES_NORM_THRESHOLD = 1e-5;

        C_DEGREE_INC = 0;
        VERIFICATION_C_DEGREE_INC = 0;
        
        VERIFICATION_PHY_DEGREE_INC = 0;
        VERIFICATION_LAMBDA_DEGREE_INC = 0;
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

