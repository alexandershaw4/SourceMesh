% Script to generate a subject specific brain surface and add overlays,
% networks, nodes and labels from the AAL90 atlas.
%
%
% Dependencies: fieldtrip & spm
% AS17

% Installation: addpath to wherever the toolbox is:
addpath(genpath('~/Downloads/MeshAAL-master'));


% Generate mesh: load, align, segment, isosurface, smooth, return gifti
g = Vol2SurfAS('my-coreg-ctf-mri.mri','mri',0.15);

% plot the mesh brain
%----------------------------------------------------------------------
atemplate('gifti',g);

% plot & save gifti file
atemplate('gifti',g,'write','MyGifti');

% for an overlay, you need a 1x90 vector where each of the 90 elements
% correspond to the 90 AAL regions.
%----------------------------------------------------------------------

% e.g. a pretend vector of t-values
overlay = randi([-2 4],90,1);

% to plot:
atemplate('gifti',g,'overlay',overlay); bigimg;

% to plot and save both gifti-mesh and gifti-overlay files:
atemplate('gifti',g,'overlay',overlay,'write', 'MyGiftiFile');

% to plot with labels
atemplate('gifti',g,'overlay',overlay,'labels'); bigimg;

% to use a template rather than subject specific surface:
atemplate('overlay',overlay,'labels');


% for a set of nodes, with only the specified nodes labelled
%----------------------------------------------------------------------
N = randi([0 1],90,1); % a binary 1x90 vector, with 1s for nodes to show
                       % do load('labels') for list

atemplate('gifti',g,'nodes',N,'labels');


% for a set of edges and the connected nodes
%----------------------------------------------------------------------
A = randi([0 10],90,90); % a 90x90 connectivity matrix

atemplate('gifti',g,'network',A,'labels');


% in development:
% add a set of tracks, as loaded with along-tract-stats toolbox
%----------------------------------------------------------------------
atemplate('gifti',g,'tracks',tracks,header);