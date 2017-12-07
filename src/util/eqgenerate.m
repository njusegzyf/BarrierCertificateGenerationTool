function [ Aeq, beq ] = eqgenerate( indvars, decvars, polyexpr)
%eqgenerate Summary of this function goes here
%   indvars: independent variables, such as: x1 x2 ...
%   decvars: decision variables, such as: p, c_alpha_beta,...
%   both independent variables and decision variables should be 1-dimension vector
%   polyexpr: polynomial expression representing the original constraint
%
% ����eqgenerate�������Ƕ԰������������Ķ���ʽ���ʽ��ȡ���߱���Լ�������Ծ��߱���Լ���������Ax = b�еĲ���A��b��
% ������indvars����ʾ���������ķ��ű���������
% decvars����ʾ���߱����ķ��ű���������
% polyexpr����ʾ����ʽ���ʽ�ķ��ű������ʽ
% ����ֵ��Aeq����ʾA�ľ���
% beq����ʾb����������
%

% in lp3
% charindvars = converttochar(indvars);
% chardecvars = converttochar(decvars);

% changed in lp4
import lp4util.symbolArrayToString
charindvars = symbolArrayToString(indvars);

% ����mupad��expand����������ʽչ��
expr = feval(symengine,'expand',polyexpr); % expand the polynomial

% ����mupad��collect����������ʽ���ݶ�����������
expr = feval(symengine,'collect',expr,charindvars); % treat independent variables as the actual polynomial variables and reoganize the polynomial

% ����mupad��poly2list����������ʽ�ֽ�Ϊ���ݶ����������ֵĸ�������ʽԪ��
coefindmon = feval(symengine,'poly2list',expr,charindvars);

% ��ȡÿ������ʽԪ�ص�ϵ����ϵ���а������߱��������ھ��߱�����Լ������Щϵ����Ϊ0�����Թ��е���ʽ������Լ��
for k = 1:1:length(coefindmon)
    dummyvar = reshape(coefindmon(k),2,1);
    coefind(k) = dummyvar(1);                  % coefficient of independent vars in constraints
    %mon(k,:)= double(dummyvar(2));         % the monimial degree, no use here
end

Aeq = zeros(length(coefind), length(decvars));
beq = zeros(length(coefind), 1);

% ��ȡϵ��Լ���ж��ھ��߱�����ϵ��
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
    %         if (mon(s,:) == zeros(1,length(mon(s,:)))) % �������Ϊ0�����ǳ�������뵽beq��
    %             beq(k,1) = -coefdec(s);           % the constant term
    %         else                                       % �������Ϊ��Ϊ0������뵽Aeq��
    %             Aeq_k = Aeq_k + coefdec(s) * mon(s,:);
    %             Aeq(k,:) = Aeq_k;
    %         end
    %     end
    
    [A,b]=equationsToMatrix(coefind(k),decvars);
    Aeq(k,:) = A;
    beq(k,1) = b;
end

end

