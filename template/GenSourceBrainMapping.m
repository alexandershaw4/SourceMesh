function MeshVertexParcelID = GenSourceBrainMapping(m,v,vi)
% Generate a pre-reigstered mapping from atlas vertices to brain surface
%
% m = mesh structure of brain surface (vertices and afces)
% v = source/atlas positions (nx3)
% vi = vector indexing which parcel/region each element of v belongs to
%

for i = 1:length(m.vertices)
    this  = cdist(m.vertices(i,:),v);
    [~,I] = min(this);
    MeshVertexParcelID(i) = vi(I);
end