classdef Lp4Config
    %LP4CONFIG Stores the configurations.
    
    properties (Constant)
        IS_DEBUG = true;
        IS_VERBOSE = false;
        
        % ���Ϊ true����֤ʱ�����˳��Դ�������� lambda �� phy�������Դ�������� phy �� lambda
        IS_VERIFY_WITH_PHY = true;
        
        % ���� rou <= ROU_THRESHOLD �ж��Ƿ��н�
        ROU_THRESHOLD = 1e-7
        
        % �����ʱ��ÿ��ϵ�������� 2 ��������ֵ
        RES_NORM_THRESHOLD = 1e-5;
        
        % �Ƿ��ӡʧ�ܵ���֤����ϸ��Ϣ
        IS_PRINT_FAILED_VERIFICATION_INFO = false;
        
        % ��һ�����ʱ��ÿ����ʽ�ұ���Ĵ������� = �����Ĵ��� + C_DEGREE_INC
        C_DEGREE_INC = 0;
        % �ڶ�����֤ʱ��ÿ����ʽ�ұ���Ĵ������� = �����Ĵ��� + VERIFICATION_C_DEGREE_INC
        VERIFICATION_C_DEGREE_INC = 0;
        
        % �ڶ�����֤ʱ���������� phy������������� = ��һ�����ʱ phy �Ĵ������� + VERIFICATION_PHY_DEGREE_INC
        VERIFICATION_PHY_DEGREE_INC = 0;
        % �ڶ�����֤ʱ���������� lambda������������� = ��һ�����ʱ lambda���Ĵ������� + VERIFICATION_PHY_DEGREE_INC
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

