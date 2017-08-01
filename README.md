# MeshAAL

Plot AAL overlays and networks (nodes & edges) on template brain natively in matlab.

```
template('labels','network',A);    % plot nodes and edges (A=90x90 double), add AAL labels 
template('labels','overlay',L);    % plot overlay (L=1x90 double), add labels
template('overlay',L,'network',A); % plot overlay and network, no labels
```

# Individual functions

Plot AAL nodes and edges on template brain natively in matlab

```
templatemesh(A);
```

where A is a 90x90 connectivity matrix of the 90 AAL nodes.


![alt text](example.gif)


Also project a 1x90 vector of node values as a mesh overlay.

```
[M,S] = templateoverlay(L) % first call
templateoverlay(L'*M,S)    % second call
```

* L is a vector of length 90 corresponding to the AAL90 atlas
* L is mapped to the size of the template by finding the n-closest points and linearly interpreting to generate a smooth surface
* Returns matrix M of weights, so that it needn't be recomputed.

![alt text](NodePowOnSurface.gif)
