function quickalign(v0,v1)

M = diag(mean([min(v0);max(v0)]...
    ./ [min(v1);max(v1)],1));

v1 = v1*M;

% singulars
[u,s,v] = svd(v0 - repmat(spherefit(v0),[size(v0,1) 1]),'econ');
[ux,sx,vx] = svd(v1 - repmat(spherefit(v1),[size(v1,1) 1]),'econ');

quickscatter3(v0,v1*(v./vx))

end