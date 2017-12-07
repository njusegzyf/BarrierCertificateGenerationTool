classdef Lp
    %LP 
    
    properties
        indvars % ����1�еĶ�������������x = [x1, x2, ..., xn]������Ϊ���Ż��������ɵ�������
        degree % ����1�еĴ������ʽ�����յĴ���������Ϊ������
        phy % ����1�еĴ������ʽ�����գ���lp.degree��ʼ�����Զ���ֵΪϵ��Ϊ���߱����Ķ���ʽ
        r %����1�еĦˡ�
        eps % ����1�е�?1��?1���ɵ�����
        f % ����1�е�f
        decvars % ����1�еľ��߱���������P, C��,��, C��,��, Cu,v
        expr
    end
    
    methods
      
      function lp = lpr(lp, r)
        lp.r = r;
      end
      
      function lp = lpf(lp, f)
        lp.f = f;
      end
      
      function lp = lpeps(lp, eps)
        lp.eps = eps;
      end
    end
    
end

