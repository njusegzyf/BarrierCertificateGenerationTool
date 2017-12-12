classdef Lp4Config
    %LP4CONFIG Stores the configurations.
    
    properties (Constant)
        IS_DEBUG = true;
        IS_VERBOSE = false;
        
        % ���Ϊ true����֤ʱ�����˳��Դ�������� lambda �� phy�������Դ�������� phy �� lambda
        IS_VERIFY_WITH_PHY = true;
        
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
        
        function displayDelimiterLine()
            disp('--------------------------------------------------------------');
        end
    end
end

