function [decvars, expression, indexVectors] = getConstraintExpressionAsVectors(degree, g, c)
% 基于 constraintExpression 修改，多返回 indexVector 数组 （所有单项式中各个约束的次数）

% expr表示要生成的多项式表达式
global expr;
expr = 0;

% indexVector表示每个单项式中各个约束的次数
global indexVector;
indexVector = -ones(1,2*length(g));

% number表示需要用到的系数变量的个数
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
% 通过递归的方式生成表达式
% 主要思想是遍历所有的次数可能性，即遍历indexVector，然后计算出indexVector表示的单项式次数，
% 如果小于给定次数degree，则加入到expr中

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
        % 优化步骤，即如果当前次数已经大于degree，则直接砍掉树的这个分支，减小开销
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