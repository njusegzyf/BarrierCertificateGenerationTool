function [decvars, expression] = constraintExpression(degree,g,c)
% ���ɴ���������degree����Լ������g��ϳɵĶ���ʽ
% ������degree������ʽ������������
% g����ʾԼ�����ϵķ��ű�������ʽ������
% c����ʾ��ʼϵ�������ķ��ű���������
%
% ����ֵ��decvars����ʾ�����ĳ�ʼϵ�������ķ��ű���������
% expression����g��Լ����ϳɵı��ʽ

% expr��ʾҪ���ɵĶ���ʽ���ʽ
global expr;
expr = 0;

% indexVector��ʾÿ������ʽ�и���Լ���Ĵ���
global indexVector;
indexVector = -ones(1,2*length(g));

% number��ʾ��Ҫ�õ���ϵ�������ĸ���
global number; % number is the index of actual variable
number = 0;

addMonomial(2*length(g),degree,g,c,0);

expression = expr;
decvars = c(1:number);

clear global expr indexVector number;

end

%% subfunction
function addMonomial(i, degree, g, c, monomialDegree)
% ͨ���ݹ�ķ�ʽ���ɱ��ʽ
% ��Ҫ˼���Ǳ������еĴ��������ԣ�������indexVector��Ȼ������indexVector��ʾ�ĵ���ʽ������
% ���С�ڸ�������degree������뵽expr��

global expr;
global indexVector;
global number;

if (i == 0)
    number = number + 1;
    mono = c(number);
    for j = 1:1:length(g)
        mono = mono * g(j)^indexVector(j);
    end
    for j = 1:1:length(g)
        mono = mono * (1-g(j))^indexVector(length(g)+j);
    end
    expr = expr + mono;
else
    for j = 0:1:degree
        % �Ż����裬�������ǰ�����Ѿ�����degree����ֱ�ӿ������������֧����С����
        indexVector(i) = j;
        %         indexVector
        %         monomialDegree
        %         i
        % optimization part, it the monimal degree is larger than degree, cut this bruch
        if (i>length(g))
            monoDegree = monomialDegree + feval(symengine, 'degree', 1-g(i-length(g))) * indexVector(i);
        end
        if (i<=length(g))
            monoDegree = monomialDegree + feval(symengine, 'degree', g(i)) * indexVector(i);
        end
        if (monoDegree>degree)
            return;
        end
        addMonomial(i-1,degree,g,c,monoDegree);
    end
end

end