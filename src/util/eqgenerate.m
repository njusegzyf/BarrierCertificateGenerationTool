function [ Aeq, beq ] = eqgenerate( indvars, decvars, polyexpr)
%eqgenerate Summary of this function goes here
%   indvars: independent variables, such as: x1 x2 ...
%   decvars: decision variables, such as: p, c_alpha_beta,...
%   both independent variables and decision variables should be 1-dimension vector
%   polyexpr: polynomial expression representing the original constraint
%
% 函数eqgenerate的作用是对包含独立变量的多项式表达式提取决策变量约束，并对决策变量约束获得形如Ax = b中的参数A和b。
% 参数：indvars：表示独立变量的符号变量向量；
% decvars：表示决策变量的符号变量向量；
% polyexpr：表示多项式表达式的符号变量表达式
% 返回值：Aeq：表示A的矩阵；
% beq：表示b的列向量。
%

% in lp3
% charindvars = converttochar(indvars);
% chardecvars = converttochar(decvars);

% changed in lp4
import lp4util.symbolArrayToString
charindvars = symbolArrayToString(indvars);

% 利用mupad中expand函数将多项式展开
expr = feval(symengine,'expand',polyexpr); % expand the polynomial

% 利用mupad中collect函数将多项式根据独立变量整理
expr = feval(symengine,'collect',expr,charindvars); % treat independent variables as the actual polynomial variables and reoganize the polynomial

% 利用mupad中poly2list函数将多项式分解为根据独立变量划分的各个单项式元素
coefindmon = feval(symengine,'poly2list',expr,charindvars);

% 获取每个单项式元素的系数，系数中包含决策变量。对于决策变量的约束是这些系数都为0。所以共有单项式个数个约束
for k = 1:1:length(coefindmon)
    dummyvar = reshape(coefindmon(k),2,1);
    coefind(k) = dummyvar(1);                  % coefficient of independent vars in constraints
    %mon(k,:)= double(dummyvar(2));         % the monimial degree, no use here
end

Aeq = zeros(length(coefind), length(decvars));
beq = zeros(length(coefind), 1);

% 提取系数约束中对于决策变量的系数
for  k = 1:1:length(coefind)
    import lp4.Lp4Config
    if Lp4Config.isVerbose()
        disp(['k:',num2str(k)]);
    end
    
    %     coefdecmon = feval(symengine,'poly2list',coefind(k),chardecvars);
    %     Aeq_k=zeros(1,length(decvars));
    %     for s = 1:1:length(coefdecmon)
    %         disp(['s:',num2str(s)]);
    %         dummyvar = reshape(coefdecmon(s),2,1);
    %         coefdec(s) = dummyvar(1);                  % coefficient of decision vars in constraints
    %         mon(s,:)= double(dummyvar(2));         % the monimial degree
    %         if (mon(s,:) == zeros(1,length(mon(s,:)))) % 如果次数为0，则是常数项，加入到beq中
    %             beq(k,1) = -coefdec(s);           % the constant term
    %         else                                       % 如果次数为不为0，则加入到Aeq中
    %             Aeq_k = Aeq_k + coefdec(s) * mon(s,:);
    %             Aeq(k,:) = Aeq_k;
    %         end
    %     end
    
    [A,b]=equationsToMatrix(coefind(k),decvars);
    Aeq(k,:) = A;
    beq(k,1) = b;
end

end

