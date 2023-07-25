function p = pca_on_mesh_overlay(D)
% compute pca on output of a sourmesh/atemplate model
% input struct returned by atemplate call, e.g. D = atemplate(' ... ')

A = spm_mesh_adjacency(D.mesh);
Y = D.overlay.data;
M = Y.*A.*Y';

[u,s,v] = svds(M);

d = max(length(s),6);

for i = 1:d
    p(i,:) = u(:,i)'*M;
end

end