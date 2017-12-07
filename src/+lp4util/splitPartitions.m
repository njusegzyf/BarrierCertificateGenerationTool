function res = splitPartitions(partitions)

% TODO

partitionsLength = length(partitions);
splitedPartitions = cell(1, partitionsLength);

for i = 1 : 1 : partitionsLength
    splitedPartitions{i} = partitions(i).splitPartition();
end

combinationCount = 1;
for i = 1 : 1 : partitionsLength
    combinationCount = combinationCount * length(splitedPartitions{i});
end

combinationCount

regionInterval = length(splitedPartitions{0});
regionLength = 1;

for i = 1 : 1 : partitionsLength
    splitedPartitionsI = splitedPartitions{i};
    splitedPartitionsILength = length(splitedPartitionsI);
    for j = 1 : 1 : splitedPartitionsILength
      splitedPartitionsICandidate = splitedPartitionsI(j);
    end
end

end
