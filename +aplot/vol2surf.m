function [y,data] = vol2surf(vol,data,wb)
% Convert a FUNCTIONAL volume to surface. Delauay triangulate, compute a
% set of 'source' vertices. Take liberties with downsampling and smoothing


% bounds:
S = size(vol);

% check if it's a 'full' volume!
if length(find(vol)) == prod(S)
    vol = vol - mode(vol(:));
end

% waitbar
waitbar(.15,wb,'Reading volume: Smoothing volume');

% a little smoothing
vol = smooth3(vol,'gaussian');

% New --- only exists is vol was loaded by atemplate
try
    pixdim = data.volume.hdr.dime.pixdim(2:4);
end

if ~isfield(data,'volume') || ~isfield(data.volume,'grid')
    fprintf('Using computed voxel coordinates\n');
    x = 1:size(vol,1);
    y = 1:size(vol,2);
    z = 1:size(vol,3);
elseif isfield(data.volume,'grid')
    fprintf('Using real-world voxel coordintes from nifti\n');
    x = data.volume.grid.x;
    y = data.volume.grid.y;
    z = data.volume.grid.z;
end

% waitbar
waitbar(.30,wb,'Reading volume: Indexing volume');

% find indices of tissue in old grid
[nix,niy,niz] = ind2sub(size(vol),find(vol));
[~,~,C]       = find(vol);

% waitbar
waitbar(.45,wb,'Reading volume: Compiling new vertex list');

% compile a new vertex list
fprintf('Compiling new vertex list (%d verts)\n',length(nix));
v = [x(nix); y(niy); z(niz)]';
v = double(v);
try
    v = v*diag(pixdim);
end

% apply affine if req.
if isfield(data.overlay,'affine')
    affine = data.overlay.affine;
    if length(affine) == 4
        fprintf('Applying affine transform\n');
        va = [v ones(length(v),1)]*affine;
        v  = va(:,1:3);
    end
end

% Fit this gridded-volume inside the box extremes of the mesh
B        = [min(data.mesh.vertices); max(data.mesh.vertices)];
V        = v - repmat(spherefit(v),[size(v,1),1]);
V(:,1)   = B(1,1) + ((B(2,1)-B(1,1))).*(V(:,1) - min(V(:,1)))./(max(V(:,1)) - min(V(:,1)));
V(:,2)   = B(1,2) + ((B(2,2)-B(1,2))).*(V(:,2) - min(V(:,2)))./(max(V(:,2)) - min(V(:,2)));
V(:,3)   = B(1,3) + ((B(2,3)-B(1,3))).*(V(:,3) - min(V(:,3)))./(max(V(:,3)) - min(V(:,3)));
v        = V;


% % new grid
% fprintf('Generating grid for volume data\n');
% x = linspace(B(1,1),B(2,1),S(1));
% y = linspace(B(1,2),B(2,2),S(2));
% z = linspace(B(1,3),B(2,3),S(3));
% 
% % find indiced of tissue in old grid
% [nix,niy,niz] = ind2sub(size(vol),find(vol));
% [~,~,C]       = find(vol);
% 
% % compile a new vertex list
% fprintf('Compiling new vertex list (%d verts)\n',length(nix));
% v = [x(nix); y(niy); z(niz)]';
% v = double(v);


% reduce patch
fprintf('Reducing patch density\n');

% waitbar
waitbar(.50,wb,'Reading volume: Triangulating');

nv  = length(v);
tri = delaunay(v(:,1),v(:,2),v(:,3));
fv  = struct('faces',tri,'vertices',v);
count  = 0;

% waitbar
waitbar(.60,wb,'Reading volume: Smoothing');

% smooth overlay at triangulated points first
Cbound = [min(C) max(C)];
C      = spm_mesh_smooth(fv,double(C),4);
C      = Cbound(1) + (Cbound(2)-Cbound(1)).*(C - min(C))./(max(C) - min(C));

% waitbar
waitbar(.70,wb,'Reading volume: Reducing patch density');

while nv > 10000
   fv  = reducepatch(fv, 0.5);
   nv  = length(fv.vertices);
   count = count + 1;
end

% print
fprintf('Patch reduction finished\n');
fprintf('Using nifti volume as sourcemodel and overlay!\n');
fprintf('New sourcemodel has %d vertices\n',nv);

% waitbar
waitbar(.90,wb,'Reading volume: Computing colours for reduced patch');

% find the indices of the retained vertexes only
fprintf('Retrieving vertex colours\n');
Ci = aplot.compute_reduced_indices(v, fv.vertices);

% waitbar
waitbar(1,wb,'Reading volume: eComplete');
close(wb);


% Update sourcemodel and ovelray data
v                    = fv.vertices;
data.sourcemodel.pos = v;
y                    = C(Ci);

end