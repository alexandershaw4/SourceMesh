function [edge,node] = rw_edgenode(name,rw,varargin)
% Read / write (AAL) edge and node files from a matrix
%
% usage:
%   read:  [edge,node] = rw_edgenode('name');
%   write: rw_edgenode('savename','w',A)
%
% where A is connectivity matrix.
% 
% To only write a subset of nodes / edges, see write_plot_node.m
% AS

if nargin < 2; rw = 'r' ; end

switch rw 
    case 'r';
    [fp,f,e] = fileparts(name);
    if isempty(fp)
        node    = dlmread([f '.node']);
        edge    = dlmread([f '.edge']);
    else
        node    = dlmread([fp '/' f '.node']);
        edge    = dlmread([fp '/' f '.edge']);
    end
    
    case 'w';
        A = varargin{1};
        A(isnan(A)) = 0;
        conmat2nodes(A,name);
end