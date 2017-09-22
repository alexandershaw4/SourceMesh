# MeshAAL

Plot functional overlays and networks (nodes & edges) from the AAL90 on template or subject cortical surface, natively in matlab. Also includes labelling and automated routine for conversion of CTF .mri file to cortical mesh.

Includes some compiled cpp code as (linux) mex files for speed up.

Designed for use with data from AAL90 beamformer MEG data.

![rotatebrains](DualRotate.gif)

# Installation 
addpath to wherever the toolbox is:
```
addpath(genpath('~/Downloads/MeshAAL-master'));
```

Dependencies: fieldtrip & spm

# Usage examples
```
  atemplate('labels');        % template mesh with AAL labels
  atemplate('overlay',L);     % template mesh with overlay
  atemplate('network',A);     % template mesh with nodes & edges
  
  atemplate('overlay',L,'network',A,'labels'); % overlay, network & labels
  atemplate('tracks',tracks,header); % plot tracks loaded with trk_read
  atemplate('gifti',g);              % use supplied gifti surface / mesh 

  atemplate('gifti'  ,g,'write',name);  % write mesh gifti
  atemplate('overlay',L,'write',name);  % write mesh & overlay giftis

```

![alt text](ExampleTracksNodesLabels.gif)



# Generate mesh from (coregistered) CTF MRI file
Auto load, align, segment, isosurface, smooth, return gifti.

Uses fieldtrip functions, isosurface & some included funcs.
```
g = Vol2SurfAS('my-coreg-ctf-mri.mri','ctf',0.15); % file, format & smooth fraction
```

Plot the mesh brain
```
atemplate('gifti',g);
```

Plot & Save (as .gii gifti file)
```
atemplate('gifti',g,'write','MyGifti');
```

![alt text](ExampleMeshRotate.gif)

# Overlay
For an overlay, you need a 1x90 vector where each of the 90 elements correspond to the 90 AAL regions.

e.g. a pretend vector of t-values:
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

Both overlay & network with labels:

![alt text](Examples/Ex_OverlayWithNetworkAndLabels.gif)

# Video
For a video from the frames in double matrix V, when V(90,n) and n is the number of frames / images that make the video.

```
m = randi([0 4],90,100);
savename = 'MyVideo';
times    = 1:1:100;

atemplate('gifti',g,'video',m,savename,times);   
``` 
![alt text](VideoExample.gif)

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

# Other tips & hints
```
atemplate('network',edges,'nosurf')   % nosurf option - no brain / template plotted.
atemplate('network',edges,'fighnd',h) % fighnd option - pass a figure handle (e.g. for subplot) 
atemplate('overlay',L    ,'nocolbar') % don't add a colorbar (e.g. overlay or network)

slice3() % copy figure to new figure with 3 subplots showing left, topography and right views
bigimg;  % make figure window bigger and orientate brain (0,0) in Az/El.
```

# in development:
add a set of tracks, as loaded with along-tract-stats toolbox
```
atemplate('gifti',g,'tracks',tracks,header);
```


