classdef LinearProgram4VerificationBase  < lp4.Lp4AndHlpVerificationBase
    
    properties
        f % ����1�е�f
        
        phy % ����1�еĶ���ʽ��
        lambda % ����1�еĶ���ʽ��
        
        c1Length % ��¼Լ�� 1-3 �ľ��߱��� C �ĳ���
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
