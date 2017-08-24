# MeshAAL

Plot AAL overlays and networks (nodes & edges) on template brain natively in matlab.

Includes some compiled cpp code as (linux) mex for speed up.

Usages:
```
  atemplate('labels');        % template mesh with AAL labels
  atemplate('overlay',L);     % template mesh with overlay
  atemplate('network',A);     % template mesh with nodes & edges
  atemplate('overlay',L,'network',A,'labels'); % overlay, network & labels
  atemplate('tracks',tracks,header); % plot tracks loaded with trk_read
  atemplate('gifti',g);       % use supplied gifti surface / mesh 

  atemplate('gifti'  ,g,'write',name);  % write mesh gifti
  atemplate('overlay',L,'write',name);  % write mesh & overlay giftis

^Function Vol2SurfAS included
```

![alt text](ExampleTracksNodesLabels.gif)



# Installation 
addpath to wherever the toolbox is:
```
addpath(genpath('~/Downloads/MeshAAL-master'));
```

Dependencies: fieldtrip & spm


# Generate mesh
Load, align, segment, isosurface, smooth, return gifti.

Uses fieldtrip functions, isosurface & some included funcs.
```
g = Vol2SurfAS('my-coreg-ctf-mri.mri','ctf',0.15);
```

Plot the mesh brain
```
atemplate('gifti',g);
```

Plot & Save (as .gii gifti file)
```
atemplate('gifti',g,'write','MyGifti');
```

# Overlay
For an overlay, you need a 1x90 vector where each of the 90 elements correspond to the 90 AAL regions.

e.g. a pretend vector of t-values
```
overlay = randi([-2 4],90,1);
atemplate('gifti',g,'overlay',overlay); bigimg;
```

To plot and save both gifti-mesh and gifti-overlay files:
```
atemplate('gifti',g,'overlay',overlay,'write', 'MyGiftiFile');
```

To plot with labels
```
atemplate('gifti',g,'overlay',overlay,'labels'); bigimg;
```

To use a template rather than subject specific surface:
```
atemplate('overlay',overlay,'labels'); % (omit 'gifti' argument)
```

![alt text](NodePowOnSurface.gif)


# Nodes
For a set of nodes, with only the specified nodes labelled.

A binary 1x90 vector, with 1s for nodes to show.

```
load('labels')          % list of the 90 AAL labels
N = randi([0 1],90,1);  % (1,90) logical / binary list
atemplate('gifti',g,'nodes',N,'labels');
```

# Network: Edges & Nodes
For a set of edges and the connected nodes.

```
A = randi([0 10],90,90); % a 90x90 connectivity matrix
atemplate('gifti',g,'network',A,'labels');
```

![alt text](example.gif)


# in development:
add a set of tracks, as loaded with along-tract-stats toolbox
```
atemplate('gifti',g,'tracks',tracks,header);
```


