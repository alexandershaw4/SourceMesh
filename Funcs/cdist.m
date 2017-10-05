function D = cdist(mv,v) %#codegen
% find distance between the vector v and the vectors of matrix mv
%
% AS17

D = sum( ( mv - repmat( v, size(mv,1), 1 )).^2, 2);