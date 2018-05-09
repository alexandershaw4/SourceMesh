function [verts,cols,bvol] = avolreg(volume,functional)
% avolreg: register values of functional 3D volume to structural 3D volume,
% return vertices and corresponding colours as well as binary volume with
% indices of vertces in
%
%


% make them double
if ~strcmp(class(volume),'double')
    volume = double(volume);
end
if ~strcmp(class(functional),'double')
    functional = double(functional);
end


% image dimentions
dims = size(volume);
init = floor(dims/2);
Q    = @squeeze;

dims    = size(volume);
func = imresize3(functional,[dims]);

% make a grid in arbitrary units 1-100
sz  = [1 100];
xg  = linspace(sz(1),sz(2),dims(1));
yg  = linspace(sz(1),sz(2),dims(2));
zg  = linspace(sz(1),sz(2),dims(3));

nit = prod(dims);
nit = round([0:.1:1]*nit);
nit = nit(2:end);

v     = []; % list of vertices to populate
c     = []; % corresponding functional data value to this vertex
n     = 0;  % counter
bvol  = zeros(dims); % empty volume with vertex indices for reversal

for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(3)
            n = n + 1;
            if ismember(n,nit)
                fprintf('%d%% done\n',10*find(ismember(nit,n)));
            end

            curr = volume(i,j,k);
            
            if any(curr)
                v = [v ; [xg(i) yg(j) zg(k)] ];
                c = [c ; func(i,j,k)         ];
                bvol(i,j,k) = length(v);
            end
        end
    end
end

fprintf('New volume has %d vertices / voxels\n',length(v));
verts = v;
cols  = c;






% if overlay
%     
%     % box 1
%     axes(h(1)); 
%     this = fliplr( Q(regfunc(init(1),:,:)))';
%     ol   = fliplr( Q(volume (init(1),:,:)))';
%     
%     im(1) = imshow(this,'colormap',bluewhitered);
%     set(im(1), 'AlphaData', ol)
%     
%     % box 2
%     axes(h(2)); 
%     this = fliplr( Q(regfunc(:,init(2),:)))';
%     ol   = fliplr( Q(volume (:,init(2),:)))';
%     
%     im(2) = imshow(this,'colormap',bluewhitered);
%     set(im(2), 'AlphaData', ol)
%     
%     % box 3
%     axes(h(3)); 
%     this = fliplr( Q(regfunc(:,:,init(3))))';
%     ol   = fliplr( Q(volume (:,:,init(3))))';
%     
%     im(3) = imshow(this,'colormap',bluewhitered);
%     set(im(3), 'AlphaData', ol)
%     
%     % box 4: 
%     
% end
% 
% newvol = [];
% %FV = isosurface(X,Y,Z,V)
% 
% end

