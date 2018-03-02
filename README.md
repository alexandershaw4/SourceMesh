# SourceMesh

Plot MEG/EEG functional overlays, networks (nodes & edges) and more on template brain in matlab.

Provide data along with it's sourcemodel coordinates (from, e.g. fieldtrip), or use the default AAL90.

Essentially it fits the sourcemodel co-ordinates (or atlas source coords) to an MNI-space gifti surface (or template brain mesh) and allows plotting of functional overlays, nodes, labels, networks (opt. with labels) and overlay-videos on this mesh. 

It can also add atlas labels to non-atlas sourcemodel data using ICP to align the cloud points. 
Can export (after alignment / interpolation) both structural and functional gifti objects.

It can also read and plot .node and .edge files.

Includes some compiled cpp code as (linux) mex for speed up for Mac and Linux. If these cause problems, delete them.

Usages/Help:

```
%
%  MESHES:
%--------------------------------------------------------------------------
%
%  atemplate()               plot a template mesh
%  atemplate('gifti',mesh)   plot a supplied (gifti) mesh
%  atemplate('gifti',mesh, 'write',name);  plot mesh & write out gifti
%  
%
%  OVERLAYS:
%--------------------------------------------------------------------------
%
%  atemplate('overlay',L);   plot template mesh with overlay from AAL90. L is [90x1]
%
%  atemplate('sourcemodel',sormod,'overlay',L)  plot template with overlay
%  values L at sourcemodel values sormod, interpolated on surface.
%
%  atemplate('gifti',mesh,'sourcemodel',sormod,'overlay',L)  plot the supplied 
%  gifti mesh with overlay values L at sourcemodel locations sormod interpolated 
%  on surface. Sormod is n-by-3, L is n-by-1.
%
%  atemplate('gifti',mesh,'sourcemodel',sormod,'overlay',L,'write','MYGifti')  
%  - This does the plot as above but writes out TWO gifti files:
%    1. MYGifti.gii is the gifti mesh 
%    2. MYGiftiOverlay.gii is the corresponding overlay data
%
%
%  **Note on sourcemodel option: If sourcemodel from Fieldtrip, swap x & y
%  by doing sm = [sourcemod(:,2),sourcemod(:,1),sourcemod(:,3)];
%
%
%  VIDEO OVERLAY:
%--------------------------------------------------------------------------
%
%  atemplate('gifti',g,'sourcemodel',sormod,'video',m,'name',times); where
%  - g      = the gifti surface to plot
%  - sormod = sourcemodel vertices
%  - m      = overlay values [vertices * ntimes] 
%  - name   = video savename
%  - times  = vector of titles (time values?)
%
%
%  NETWORKS:
%--------------------------------------------------------------------------
%
%  atemplate('network',A);    plot template mesh with 90x90 AAL network, A.
%
%  atemplate('sourcemodel',sormod,'network',A);  plot network A  at
%  sourcemodel locations in sormod. sormod is n-by-3, netowrk is n-by-n.
%
%  atemplate('sourcemodel',sormod,'network',A,'write','savename'); 
%   - as above but writes out .node and .edge files for the network, and
%   the gifti mesh file.
%
%
%  OTHER
%--------------------------------------------------------------------------
%
%  atemplate('labels');         plot node labels (AAL90) 
%
%  atemplate('labels', all_roi_tissueindex, labels); where all_roi_tissue 
%  is a 1-by-num-vertices vector containing indices of the roi this vertex
%  belongs to, and 'labels' contains the labels for each roi. The text
%  labels are added at the centre of the ROI.
%  
%  Labels notes:
%     - If plotting a network, only edge-connected nodes are labelled.
%     - If plotting a set of nodes (below), only those are labelled.
%     - Otherwise, all ROIs/node labels are added!
%
%  atemplate('nodes', N);             Plot dots at node==1, i.e. N=[90,1]
%  atemplate('tracks',tracks,header); plot tracks loaded with trk_read
%
%  Note: any combination of the inputs should be possible.
%
%
%
%
%  AN EXAMPLE NETWORK: from 5061 vertex sourcemodel with AAL90 labels
%--------------------------------------------------------------------------
%
% load New_AALROI_6mm.mat       % load ft source model, labels and roi_inds
%
% net  = randi([0 1],5061,5061);   % generate a network for this sourmod
% pos  = template_sourcemodel.pos; % get sourcemodel vertices
% labs = AAL_Labels;               % roi labels
% rois = all_roi_tissueindex;      % roi vertex indices
%
% atemplate('sourcemodel',pos,'network',net,'labels',rois,labs);
%
%
%
%
%
%  Cortical mesh from mri
%  ---------------------------
%  If the gifti option is included, input (g) may be the filename of a 
%  coregistered ctf .mri file. This will call Vol2SurfAS which uses 
%  fieldtrip & isosurface to normalise, align, segment and extract a
%  cortical surface. This is then centred, smoothed and converted to a
%  gifti object.
%
%  See also: slice3() slice2()
%
% ^trk_read requires along-tract-stats toolbox
%
% AS17

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

![alt text](ExampleMeshRotate.gif)

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


# in development:
add a set of tracks, as loaded with along-tract-stats toolbox
```
atemplate('gifti',g,'tracks',tracks,header);
```


