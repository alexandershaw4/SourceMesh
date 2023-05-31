function [f,nv,Ci] = points_to_alphashape(v)
% given a nx3 vertex list, this function
%
% 1) triangulates faces using alphaShape which outperforms delaunay
% 2) also does minor patch reduction on the fly
% 3) computed indices of vertices retained from original list (Ci)
%
% [f,nv,Ci] = points_to_alphashape(v)
%
% AS

tr = alphaShape(v);
[f, nv] = boundaryFacets(tr);
Ci = compute_reduced_indices(v,nv);