function res = evaluatePointWithBarrier(barrier, symbols, point)

res = subs(barrier, symbols, point);

end

