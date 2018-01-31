function res = reshapeToVector(input)

if size(input, 1) ~= 1
    res = reshape(input, 1, size(input, 1) * size(input, 2));
else
    res = input;
end

end
