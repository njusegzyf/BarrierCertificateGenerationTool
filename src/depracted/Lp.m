classdef Lp
    %LP 
    
    properties
        indvars % 问题1中的独立变量，例如x = [x1, x2, ..., xn]，类型为符号化变量构成的行向量
        degree % 问题1中的带求多项式函数φ的次数，类型为正整数
        phy % 问题1中的带求多项式函数φ，当lp.degree初始化后，自动赋值为系数为决策变量的多项式
        r %问题1中的λ。
        eps % 问题1中的?1和?1构成的向量
        f % 问题1中的f
        decvars % 问题1中的决策变量，包括P, Cα,β, Cγ,δ, Cu,v
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

