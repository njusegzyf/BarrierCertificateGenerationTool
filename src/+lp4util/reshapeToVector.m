function res = reshapeToVector(input)
  res = reshape(input, 1, size(input, 1) * size(input, 2));
end

