function atlas = interp_template(atlas,pos)
% The Euclidean/ICP routine for atlas registration

if length(atlas.pos) == length(pos)
    fprintf('Overlay and atlas Vectors already match!\n');
    atlas.M = eye(length(pos));
    return;
end

fprintf('Scanning points:\n');
M = zeros( length(atlas.pos), length(pos) )';
r = floor( length(atlas.pos) / length(pos) );%1;
w = 1;

dist  = cdist(pos,atlas.pos)';    
%for i = 1:length(atlas.pos)
for i = 1:length(pos)    
    if i > 1; fprintf(repmat('\b',[size(str)])); end
    str = sprintf('%d/%d',i,(length(pos)));
    fprintf(str);

    [junk,ind] = maxpoints(dist(:,i),r,'min');
    M (i,ind)  = w;
end
fprintf('\n');
atlas.M = M';

end