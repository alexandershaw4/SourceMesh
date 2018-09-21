function out = ug(mesh)
% un-gifti object: convert gifti to struct

out          = struct;
out.vertices = mesh.vertices;
out.faces    = mesh.faces;