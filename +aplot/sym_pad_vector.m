function y = sym_pad_vector(x,n)
% Symmetrically pad a vector with zeros

if length(x) ~= n
    k = n - length(x);
    k = floor(k/2);
    y = [zeros(1,k) x zeros(1,k)];
    
else y = x;
end

end