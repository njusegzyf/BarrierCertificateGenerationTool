function number = monomialNumber(variableNumber, degree)
%   monomialNumber returns the number of monomials which consists of 'variableNumber' variables and the degree is less than 'degree'
%   函数monimialNumber的作用是根据变量的个数和degree计算出由变量构成的次数不大于degree的所有单项式个数。
%   参数：variableNumber：正整数，表示变量个数；
%   degree：正整数，表示单项式的最高允许次数
%   返回值：number：单项式个数。

number = 0;
for i = 0 : degree
    % nchoosek(variableNumber+i-1, variableNumber-1) 计算 总次数为 i 时的多项式个数
    % 这个问题是排列组合中的 N个小球放M个盒子，球同，盒不同，允许空箱 问题，用 插板法 可解
    % @see http://chensmiles.blog.163.com/blog/static/121463991200962113136292/
    number = number + nchoosek(variableNumber+i-1, variableNumber-1);
end

end

