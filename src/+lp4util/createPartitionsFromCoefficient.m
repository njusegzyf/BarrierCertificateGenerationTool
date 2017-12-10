function res = createPartitionsFromCoefficient(coefficient, delta)
if size(coefficient, 2) ~= 1
    error('');
end

coefficientSize = size(coefficient, 1);
import lp4util.Partition
res = repmat(Partition(0, 0), 1, coefficientSize);

for i = 1 : coefficientSize
    midValue = coefficient(i, 1);
    res(i) =  Partition.createPartitionFromMidValue(midValue, delta);
end

end
