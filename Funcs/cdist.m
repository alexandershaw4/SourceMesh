function D = cdist(mv,v) %#codegen
% find distance between the 3D coordiates specified in vector v (1x3) and 
% the vectors of matrix mv (nx3), or between matrix v (mx3) and matrix mv
% (nx3), where v(1,3) = [x y z]
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
