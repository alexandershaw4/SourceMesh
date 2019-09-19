function mesh = read_nv(X)

if nargin < 1
    %fprintf('Using template\n');
    X = dlmread('BrainMesh_ICBM152_smoothed.nv');
else
    if ~strcmp(X(end-1:end),'nv')
        X = [X '.nv'];
    end
    X = dlmread(X);
end

NVert = X(1);
%fprintf('There are %d vertices ',NVert);
v     = X(2:NVert+1,:);
NFace = X(NVert+2);
%fprintf('and %d faces\n',NFace);
f     = X(NVert+3:NVert+3+(NFace-1),:);

mesh.vertices = v;
mesh.faces    = f;