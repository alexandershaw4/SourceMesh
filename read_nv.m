function mesh = read_nv()

X     = dlmread('BrainMesh_ICBM152_smoothed.nv');
NVert = X(1);
v     = X(2:NVert+1,:);
NFace = X(NVert+2);
f     = X(NVert+3:NVert+3+(NFace-1),:);

mesh.vertices = v;
mesh.faces    = f;