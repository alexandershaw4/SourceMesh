
% Example of how to re-compute an already registered sourcemodel projection
% on the mesh:

% Original run:
%--------------------------------------------------------------------------

% sourcemodel with 5061 voxels/vertices
load New_AALROI_6mm.mat; pos = template_sourcemodel.pos;

% some overlay values for each vertex
o = randi([-4 7],5061,1);

% make the original plot
figure,
A = atemplate('overlay',o,'sourcemodel',pos)



% using the returned data, plot some different values from the same
% sourcemodel (much quicker):
%--------------------------------------------------------------------------

new_overlay_data = o*-1;


% get stuff returned from preivous call
M    = A.overlay.smooth_weights;
indz = A.overlay.indz;

% incorporate overlay into precomputed weights matrix
for k = 1:length(new_overlay_data)
    M(k,indz(k,:)) = new_overlay_data(k) * M(k,indz(k,:));
end

% normalise by number of overlapping points at this vertex
for i = 1:size(M,2)
    y(i) = sum( M(:,i) ) / length(find(M(:,i))) ;
end

% rescale y by L limits
S  = [min(new_overlay_data(:)),max(new_overlay_data(:))];
y  = S(1) + ((S(2)-S(1))).*(y - min(y))./(max(y) - min(y));
y(isnan(y)) = 0;

% do plot:
figure, A2 = atemplate('overlay',y)
