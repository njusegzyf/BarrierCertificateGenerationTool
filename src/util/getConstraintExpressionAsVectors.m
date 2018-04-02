function [decvars, expression, indexVectors] = getConstraintExpressionAsVectors(degree, g, c)
% ���� constraintExpression �޸ģ��෵�� indexVector ���� �����е���ʽ�и���Լ���Ĵ�����

% expr��ʾҪ���ɵĶ���ʽ���ʽ
global expr;
expr = 0;

% indexVector��ʾÿ������ʽ�и���Լ���Ĵ���
global indexVector;
indexVector = -ones(1,2*length(g));

% number��ʾ��Ҫ�õ���ϵ�������ĸ���
global number; % number is the index of actual variable
number = 0;

global globalIndexVectors;
globalIndexVectors = [];

addMonomial(2*length(g),degree,g,c,0);

expression = expr;
decvars = c(1:number);

indexVectors = globalIndexVectors;

clear global expr indexVector number globalIndexVectors;

end

%% subfunction
function addMonomial(i, degree, g, c, monomialDegree)
% ͨ���ݹ�ķ�ʽ���ɱ��ʽ
% ��Ҫ˼���Ǳ������еĴ��������ԣ�������indexVector��Ȼ������indexVector��ʾ�ĵ���ʽ������
% ���С�ڸ�������degree������뵽expr��

global expr;
global indexVector;
global number;
global globalIndexVectors;

if (i == 0)
    number = number + 1;
    globalIndexVectors(number) = indexVector;
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