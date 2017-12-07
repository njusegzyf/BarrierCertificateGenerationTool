classdef SymbolicPolynomial
  %SYMBOLICPOLYNOMIAL
    
  properties
    symbolicVars
    maxDegree
    coefficientVars
    expression
  end
  
  methods

    function this = SymbolicPolynomial(symbolicVarsArg, maxDegreeArg, coefficientVarsArg, expressionArg)
      %SYMBOLICPOLYNOMIAL 构造此类的实例
      this.symbolicVars = symbolicVarsArg;
      this.maxDegree = maxDegreeArg;
      this.coefficientVars = coefficientVarsArg;
      this.expression = expressionArg;
    end

  end

  methods (Static = true)

    function res = createSymbolicPolynomial(coefficientSymbolName, symbolicVars, maxDegree)
      symbolVarsNum = length(symbolicVars);
      coefficientVars = sym(coefficientSymbolName, [1, monomialNumber(symbolVarsNum, maxDegree)]);
      expression = coefficientVars * monomials(symbolicVars, [0:maxDegree]);

      import lp4util.SymbolicPolynomial
      res = SymbolicPolynomial(symbolicVars, maxDegree, coefficientVars, expression);
    end

  end 
end

