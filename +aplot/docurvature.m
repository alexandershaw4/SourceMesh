function C = docurvature(M)
% Computes the curvature of a mesh
%
%

A = spm_mesh_adjacency(M);
A = sparse(1:size(M.vertices,1),1:size(M.vertices,1),1./sum(A,2)) * A;

C = (A-speye(size(A))) * double(M.vertices);
N = spm_mesh_normals(M);
C = sign(sum(N.*C,2)) .* sqrt(sum(C.*C,2));

end