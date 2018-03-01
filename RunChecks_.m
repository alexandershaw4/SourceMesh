function RunChecks_()

% meshes
fprintf('Testing template mesh\n');
atemplate(); close
fprintf('Done\n\n\n');

fprintf('Testing custom gifti mesh\n');
g = gifti('spm.surf.gii');
atemplate('gifti',g); close
fprintf('Done\n\n\n');

% AAL stuff

fprintf('Testing AAL90 overlays and networks etc\n');
O  = randi([0 16],90,1);
atemplate('overlay',O,'labels'); close
fprintf('Overlay done\n');
s  = 90;
N  = randi([0 9],s,s);
ncomp = 90;
while ncomp > 20
    N = N .* randi([0 1],s,s);
    [~,ncomp] = PEig90(N);
end
atemplate('network',N,'labels'); close
fprintf('Network done\n');

Nodes = randi([0 1],90,1);
atemplate('nodes',Nodes,'labels'); close
fprintf('Nodes done\n');

% CUSTOM SOURCEMOD STUFF
load('New_AALROI_6mm.mat')
pos  = template_sourcemodel.pos;
labs = AAL_Labels;
rois = all_roi_tissueindex;
O    = randi([0 16],5061,1);

fprintf('Testing custom sourcemodel with custom gifti mesh, overlay and atlas labels\n');
atemplate('gifti',g,'sourcemodel',pos,'overlay',O,'labels',rois,labs);close;
fprintf('Done\n\n\n');

s  = 5061;
N  = zeros(s,s);
r  = round(s/100):round(s/100):s(end);
r  = r(4:4:end);
N(r,r) = 1;

fprintf('Testing custom sourcemodel with custom gifti mesh, NETWORK and atlas labels\n');
atemplate('gifti',g,'sourcemodel',pos,'network',N,'labels',rois,labs);close;
fprintf('Done\n\n\n');

