function exprs = generateRamdomExprs(vars, degree, exprCount, coffStart, coffEnd, randomSeed)
% Note: `rng` accepts a int as random seed, and if we pass a double then it is convert to int.
% So if we pass a double generate by rand from 0 to 1, then the double value is always convert to 0
% and we always get an same output.
if ~isinteger(randomSeed)
    error('Error. randomSeed must be an integer, not a %s.', class(randomSeed));
end
% errorIfWrongType(randomSeed, 'int32', 'randomSeed');

exprs(1, exprCount) = vars(1);

rng(randomSeed);

exprMonomials =  monomials(vars, 0 : degree);
exprMonomialsLen = length(exprMonomials);

for i = 1 : exprCount
    exprCoffs = coffStart + (coffEnd - coffStart) * rand(1, exprMonomialsLen);
    exprs(i) = exprCoffs * exprMonomials;
end

end

