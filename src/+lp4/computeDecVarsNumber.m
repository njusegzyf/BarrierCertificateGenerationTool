function decVarNum = computeDecVarsNumber(varsN, varsD, phyD, lambdaD)
  
  decVarNum = 1; % rou
  
  phyVarNum = monomialNumber(varsN, phyD);
  decVarNum = decVarNum + phyVarNum;
  
  lambdaVarNum = monomialNumber(varsN, lambdaD); 
  decVarNum = decVarNum + lambdaVarNum;

  wVarNum = phyVarNum * lambdaVarNum;
  decVarNum = decVarNum + wVarNum;
  
  % TODO add dec vars num for three eqs
  
end

