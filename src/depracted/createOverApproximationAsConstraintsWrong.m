function constraints = createOverApproximationAsConstraintsWrong(p, pPartition, c, cPartition, w)
  %CREATEOVERAPPROXIMATION Create an over approximation for W = P * C.
  % That is, create constraints for p = [p1, p2], c = [c1, c2], w = p * c
  % which includes the left part of constraints (for `a <= 0`, returns `a`):
  % p1 <= p <= p2
  % c1 <= c <= c2
  % (p - (p1 + p2)/2)(c - (c1 + c2)/2) <= (p2 - p1) * (c2 - c1) / 4
  % (p - (p1 + p2)/2)(c - (c1 + c2)/2) >= 0 
  
  if ~isa(p, 'sym') || ~isa(c, 'sym') || ~isa(w, 'sym')
       error('p or c or w is not a symbolic variable.');
  end
  
  pMid = (pPartition.boundLow + pPartition.boundHigh) / 2;
  cMid = (cPartition.boundLow + cPartition.boundHigh) / 2;
  
%   constraintP1 = pPartition.boundLow - p;
%   constraintP2 = p - pPartition.boundHigh;
%   
%   constraintC1 = cPartition.boundLow - c;
%   constraintC2 = c - cPartition.boundHigh;
  
  constraintW1 = -w - (cMid * p) - (pMid * c)...
      - (pPartition.boundHigh - pPartition.boundLow) * (cPartition.boundHigh - cPartition.boundLow) / 4;
  constraintW2 = w + (cMid * p) + (pMid * c);
  
  % constraintP1, constraintP2, constraintC1, constraintC2, 
  constraints = [constraintW1, constraintW2];
end

function [resA, resb] = createOverApproximationAsMatrix(pPartition, cPartition)
  %CREATEOVERAPPROXIMATION Create an over approximation for W = P * C.
  
  pMid = (pPartition.boundLow + pPartition.boundHigh) / 2;
  cMid = (cPartition.boundLow + cPartition.boundHigh) / 2;
  
  % resA * [w, p, c] < b
  resA = [-1, -cMid, -pMid;
          1, cMid, pMid;
          0, 1, 0;
          0, -1, 0;
          0, 0, 1;
          0, 0, -1];
  resb = [cMid * pMid, 0, pPartition.boundHigh, -pPartition.boundLow, cPartition.boundHigh, -cPartition.boundLow];
end

function res = createExpression1(p, pBoundLow, pBoundHigh, c, cBoundLow, cBoundHigh) 
f = (p - (pBoundLow + pBoundHigh) / 2) * (c - (cBoundLow + cBoundHigh) / 2);
simpleF = expand(f);
res = subs(simpleF, p*c, "w");
end
