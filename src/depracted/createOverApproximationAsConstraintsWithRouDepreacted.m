function constraints = createOverApproximationAsConstraintsWithRouDepreacted(p, pPartition, c, cPartition, w, rou)
%CREATEOVERAPPROXIMATION Create an over approximation for W = P * C.

if ~isa(p, 'sym') || ~isa(c, 'sym') || ~isa(w, 'sym') || ~isa(rou, 'sym')
    error('p or c or w is not a symbolic variable.');
end

if pPartition.boundLow == -1 && pPartition.boundHigh == 1 && cPartition.boundLow == -1 && cPartition.boundHigh == 1
    constraints = [-w + p + c -1 - rou,...
        w - p - c - 3 - rou ,... % + 0.5,...
        w - 1 - rou,...
        -w - 1 - rou];
elseif pPartition.boundLow == -0.5 && pPartition.boundHigh == 0.5 && cPartition.boundLow == -0.5 && cPartition.boundHigh == 0.5
    constraints = [-w + 0.5*p + 0.5*c - 0.75 - rou,...
        w - 0.5*p - 0.5*c - 0.75 - rou,...
        w - 0.25 - rou,...
        -w - 0.25 - rou];
elseif pPartition.boundLow == -0.3 && pPartition.boundHigh == 0.3 && cPartition.boundLow == -0.3 && cPartition.boundHigh == 0.3
    constraints = [-w + 0.3*p + 0.3*c - 0.27 - rou,...
        w - 0.3*p - 0.3*c - 0.27 - rou,...
        w - 0.09 - rou,...
        -w - 0.09 - rou];
elseif pPartition.boundLow == -0.2 && pPartition.boundHigh == 0.2 && cPartition.boundLow == -0.2 && cPartition.boundHigh == 0.2
    constraints = [-w + 0.2*p + 0.2*c - 0.04 - rou,... %- 0.12
        w - 0.2*p - 0.2*c + 0.02 - rou,... %- 0.12
        w - 0.04 - rou,...
        -w - 0.04 - rou];
elseif pPartition.boundLow == -0.15 && pPartition.boundHigh == 0.15 && cPartition.boundLow == -0.15 && cPartition.boundHigh == 0.15
    constraints = [-w + 0.15*p + 0.15*c - 0.0675 - rou,... %- 0.12
        w - 0.1*p - 0.1*c - 0.0675 - rou,... %- 0.12
        w - 0.0225 - rou,...
        -w - 0.0225 - rou];
elseif pPartition.boundLow == -0.1 && pPartition.boundHigh == 0.1 && cPartition.boundLow == -0.1 && cPartition.boundHigh == 0.1
    constraints = [-w + 0.1*p + 0.1*c - 0.03 - rou,... %- 0.12
        w - 0.1*p - 0.1*c - 0.03 - rou,... %- 0.12
        w - 0.01 - rou,...
        -w - 0.01 - rou];
else
    error('');
end


end
