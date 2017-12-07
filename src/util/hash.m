function index = hash( indexVector, degree )
%  HASH
%  indexVector: the vector which represents the subscript index of the
%  coefficient;
%  degree: the degree of the Lyapunov function set by users
%  hash function treats the indexVector as a (degree+1)-radix number, and
%  returns the corresponding decimal number as the index. Note that if the
%  vector is in the form of [0,0,...,0], the corresponding index will be 0,
%  which is illegal. The actual return value is the corresponding decimal
%  nuber + 1.

index = 1;
for i = 1:1:length(indexVector)
    index = index + indexVector(i) * (degree + 1)^(length(indexVector) - i);
end

end

