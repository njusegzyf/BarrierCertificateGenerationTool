classdef Lp4Config
    %LP4CONFIG Stores the configurations.
    
    properties (Constant)
        IS_DEBUG = true;
        IS_VERBOSE = false;
        
        % 如果为 true，验证时，除了尝试代入求出的 lambda 求 phy，还尝试代入求出的 phy 求 lambda
        IS_VERIFY_WITH_PHY = true;
        
        % 根据 rou <= ROU_THRESHOLD 判断是否有解
        ROU_THRESHOLD = 1e-7
        
        % 检验解时，每个系数向量的 2 范数的阈值
        RES_NORM_THRESHOLD = 1e-5;
        
        % 是否打印失败的验证的详细信息
        IS_PRINT_FAILED_VERIFICATION_INFO = false;
        
        % 第一步求解时，每个公式右边项的次数限制 = 左边项的次数 + C_DEGREE_INC
        C_DEGREE_INC = 0;
        % 第二步验证时，每个公式右边项的次数限制 = 左边项的次数 + VERIFICATION_C_DEGREE_INC
        VERIFICATION_C_DEGREE_INC = 0;
        
        % 第二步验证时，如果求的是 phy，则其次数限制 = 第一步求解时 phy 的次数限制 + VERIFICATION_PHY_DEGREE_INC
        VERIFICATION_PHY_DEGREE_INC = 0;
        % 第二步验证时，如果求的是 lambda，则其次数限制 = 第一步求解时 lambda，的次数限制 + VERIFICATION_PHY_DEGREE_INC
        VERIFICATION_LAMBDA_DEGREE_INC = 0;
        
        DEFAULT_DEC_VAR_SIZE = 2048;
        
        DEFAULT_PARTITION_REPEAT_NUM = 1024;
        
        DEFAULT_EPS = 0.00001
        
        IS_SET_LINPROG_LOWERBOUND = true;
        LINPROG_LOWERBOUND = -100000000000000;
        
        IS_USE_CVX = true;
    end
    
    methods (Static)
        
        function res = processDegree(de)
            if mod(de, 2) == 0
                res = de;
            else
                res = de + 1;
            end
        end
        
        function res = isDebug()
            import lp4.Lp4Config
            res = Lp4Config.IS_DEBUG;
        end
        
        function res = isVerbose()
            import lp4.Lp4Config
            res = Lp4Config.IS_VERBOSE;
        end
        
        function displayDelimiterLine()
            disp('--------------------------------------------------------------');
        end
        
        function displaySolveRes(solveRes, resNorms)
            if ~(solveRes.hasSolution())
                disp('Verify failed.');
            else
                disp('Verify succeed, norms ;');
                disp(resNorms);
            end
        end
        
        function [lpVer, solveResVer, resNorms] = createAbsentVerificationResult()
            lpVer = 0;
            solveResVer = 0;
            resNorms = [];
        end
    end
end

