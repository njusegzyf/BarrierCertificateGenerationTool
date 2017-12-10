function constraints = createOverApproximationAsConstraints(p, pPartition, c, cPartition, w)
%createOverApproximationAsConstraints Create an over approximation for W = P * C.

if ~isa(p, 'sym') || ~isa(c, 'sym') || ~isa(w, 'sym')
    error('p or c or w is not a symbolic variable.');
end

if pPartition.boundLow == -1 && pPartition.boundHigh == 1 && cPartition.boundLow == -1 && cPartition.boundHigh == 1
    constraints = [-w + p + c -1,...
        w - p - c - 3,... % + 0.5,...
        w - 1,...
        -w - 1];
elseif pPartition.boundLow == -0.5 && pPartition.boundHigh == 0.5 && cPartition.boundLow == -0.5 && cPartition.boundHigh == 0.5
    constraints = [-w + 0.5*p + 0.5*c - 0.75,...
        w - 0.5*p - 0.5*c - 0.75,...
        w - 0.25,...
        -w - 0.25];
elseif pPartition.boundLow == -0.3 && pPartition.boundHigh == 0.3 && cPartition.boundLow == -0.3 && cPartition.boundHigh == 0.3
    constraints = [-w + 0.3*p + 0.3*c - 0.27,...
        w - 0.3*p - 0.3*c - 0.27,...
        w - 0.09,...
        -w - 0.09];
elseif pPartition.boundLow == 0 && pPartition.boundHigh == 4 && cPartition.boundLow == 0 && cPartition.boundHigh == 4
    constraints = [-w + 2*p + 2*c -4, w - 2*p - 2*c];
else
    
    constraintW1 = expand(-(p - pPartition.boundLow)*(c - pPartition.boundLow));
    constraintW1 = subs(constraintW1, p*c, w);
    
    constraintW2 = expand(-(p - pPartition.boundHigh)*(c - pPartition.boundHigh));
    constraintW2 = subs(constraintW2, p*c, w);
    
    pMid = (pPartition.boundLow + pPartition.boundHigh) / 2;
    cMid = (cPartition.boundLow + cPartition.boundHigh) / 2;
    
    constraintW3 = w - (cMid * p) - (pMid * c) + cMid * pMid...
        - (pPartition.boundHigh - pPartition.boundLow) * (cPartition.boundHigh - cPartition.boundLow) / 4;
    
    constraints = [constraintW1, constraintW2, constraintW3];
end

%   pRange = pPartition.boundHigh - pPartition.boundLow;
%   cRange = cPartition.boundHigh - cPartition.boundLow;
%   x = 4 * (p - pPartition.boundLow) / pRange;
%   y = 4 * (c - cPartition.boundLow) / cRange;
%
%   constraintW1 = expand(-x*y + 2*x + 2*y - 4);
%   constraintW1 = subs(constraintW1, p*c, w);
%
%   constraintW2 = expand(x*y - 2*x - 2*y);
%   constraintW2 = subs(constraintW2, p*c, w);
%
%   constraints = [constraintW1, constraintW2];

end
