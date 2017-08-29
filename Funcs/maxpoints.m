function [V,I] = maxpoints(x,n,varargin)
% Find max n-points in vector x, and their indices.
%
% optional third input, function string, e.g. 'min'
% hint - also use abs(x) input
%
%

if nargin < 3
    f = @max;
else
    f = eval(['@' lower(varargin{1})]);
end

for i = 1:n
   
    [V(i),I(i)] = f(x);
    x(I(i)) = nan; 
    
end