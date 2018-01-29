function exprs = generateRamdomExprs(vars, degree, exprCount, coffStart, coffEnd, randomSeed)

exprs(1, exprCount) = vars(1);

rng(randomSeed);

exprMonomials =  monomials(vars, 0 : degree);
exprMonomialsLen = length(exprMonomials);

for i = 1 : exprCount
    exprCoffs = coffStart + (coffEnd - coffStart) * rand(1, exprMonomialsLen);
    exprs(i) = exprCoffs * exprMonomials;
end

end

