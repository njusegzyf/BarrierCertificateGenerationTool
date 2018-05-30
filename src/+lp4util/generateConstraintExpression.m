function [decvars, expression] = generateConstraintExpression(maxDegree, g, cNamePrefix)
% 生成次数不大于degree的,由约束集合g组合成的多项式
% 参数：degree：多项式次数正整数；
% g：表示约束集合的符号变量多项式向量；
%
% 返回值：decvars：表示真正的初始系数变量的符号变量向量；
% expression：由g中约束组合成的表达式

gLen = length(g);
basesLen = gLen * 2;
% pre allocate memory for arrays
bases(basesLen) = g(1);
baseDegrees(basesLen) = 0;

for i = 1 : gLen
    base = g(i);
    
    % add two bases for g
    bases(2*i - 1) = base;
    bases(2*i) = 1 - base;

    % compute the degree of the bases
    baseDegree = feval(symengine, 'degree', base);
    
    baseDegrees(2*i) = baseDegree;
    baseDegrees(2*i - 1) = baseDegree;
end

basesPowers = zeros(1, basesLen);

% 预分配决策变量，最大数目为基底均为1次时使用的决策变量数目
c = sym(cNamePrefix, [1, lp4.Lp4Config.getDecVarArraySize(basesLen, maxDegree)]);

% 实际使用的决策变量数目
cNum = 1;

% init the expression with constant term
expression = c(cNum);

    function addCombinations(currentPos, currentDegree)
        
        if currentPos > basesLen
            return
        end
        
        if currentDegree <= maxDegree
            % first add the current combination if the power for current pos base is not zero
            if (basesPowers(currentPos) ~= 0)
                
                % compute the current combination
                currentCombination = 1;
                for j = 1 : currentPos
                    currentCombination = currentCombination*bases(j)^basesPowers(j);
                end
                
                % add
                cNum = cNum + 1;
                expression = expression + c(cNum) * currentCombination;
            end
            
            % if we can add more bases
            if currentDegree < maxDegree
                % add reamin bases
                % clear power for the base at next position
                basesPowers(currentPos + 1) = 0;
                addCombinations(currentPos + 1, currentDegree);
                
                % add the current base
                baseDegree = baseDegrees(currentPos);
                currentDegree = currentDegree + baseDegree;
                basesPowers(currentPos) = basesPowers(currentPos) + 1;
                addCombinations(currentPos, currentDegree);               
            end
            
        end % if current degree is already larger than degree, cut the branch
    end

% start from first base
addCombinations(1, 0);

decvars = c(1, 1:cNum);

end
