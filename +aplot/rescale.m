function y = rescale(x,S)

y = S(1) + (S(2)-S(1)) .* (x - min(x) ) / ...
    ( max(x) - min(x) );

end