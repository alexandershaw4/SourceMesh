function D = cdist(mv,v) %#codegen

D = sum( ( mv - repmat( v, size(mv,1), 1 )).^2, 2);