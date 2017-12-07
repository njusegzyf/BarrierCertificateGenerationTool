function degree = computeDegree(expr, indvars)
% 计算一个多项式的次数
% expr：含符号变量的多项式表达式
% indvars：独立变量构成的向量
% degree：表示多项式次数的正整数

% 利用mupad中poly2list函数将多项式分解为各个单项式元素。
% 每个元素由是一个1 × 2的向量，第一个元素表示该单项式的系数，第二个元素表示该单项式对于每个符号变量的次数。
% 在得到所有单项式后，计算出每个单项式的次数，即所有独立变量的次数之和，并获取最大的次数作为整个多项式的次数。

% charindvars = converttochar(indvars);
% Note: In Matlab 9.3.0.713579 (R2017b), `converttochar` returns a string
% like `(matrix([[x1, x2]]`, which does not work.

import lp4util.symbolArrayToString
charindvars = symbolArrayToString(indvars);

coefmon = feval(symengine, 'poly2list', expr, charindvars);

degree = 0;

for k = 1:1:length(coefmon)
    dummyvar = reshape(coefmon(k), 2, 1);
    
    % `mon(k,:)` defines the k-th row of the matrix mon.
    % mon(k,:) = double(dummyvar(2));
    % if (sum(mon(k,:)) > degree)
    %     degree = sum(mon(k,:));
    % end
    
    % changed to
    tempDegree = sum(double(dummyvar(2)));
    if (degree < tempDegree)
        degree = tempDegree;
    end
end

end

