function D = cdist(mv,v) %#codegen
% find distance between the vector v and the vectors of matrix mv
%
% AS17

if isvector(v)
    D = sum( ( mv - repmat( v, size(mv,1), 1 )).^2, 2);
end

if ismatrix(v)
    for i = 1:3
        dif(:,:,i) = (v(:,i)'-mv(:,i) );
    end
    D = squeeze(sum(( dif.^2),3));
end
