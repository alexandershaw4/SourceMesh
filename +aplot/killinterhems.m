function x = killinterhems(x);
% kill interhemispherics, assuming a square adjacency matrix
%
%

S  = size(x);
xb = (S(1)/2)+1:S(1);
yb = (S(2)/2)+1:S(2);
xa = 1:S(1)/2;
ya = 1:S(2)/2;

x(xa,yb) = 0;
x(xb,ya) = 0;

end