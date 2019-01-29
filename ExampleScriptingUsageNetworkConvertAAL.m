
% Plot a network on structural image, interpolated to AAL90
%==========================================================================

% from nifti T1 and network .node/.edge files:
struct = 'GreyMRI.nii';
net    = 'network.edge';

% % or from a gifti MR surface & matrices representing network (n x n) and
% % source locations (n x 3)
% struct = 'TestMesh.gii';
% net    = randi([0 1],90,90);
% funct  = (net.*net') .* randi([-7 7],90,90);
% locs   = randi([-70 70],90,3);

% parse source model locations
%--------------------------------------------------------------------------
i.A  = net;
data = aplot.sort_sourcemodel([],i);

% make struct surface
%--------------------------------------------------------------------------
i.g  = struct;
[mesh,data] = aplot.get_mesh(i,data);       


% plot structural part
%--------------------------------------------------------------------------
i.hemi      = 'both'; % if only one hemisphere plot
i.affine    = [];     % if affine supplied or flip flag
i.flip      = 0;
i.inflate   = 0;      % if inflate, pass flag
i.checkori  = 0;      % check orientation?
i.fillholes = 0;      % fill holes during hemisphere separation?
i.pmesh     = 1;      % plot the mesh or just return the info
i.write     = 0;      % write out a gifti file 
i.fname     = [];     % file name for the gifti
i.fighnd    = [];     % axes handle to plot into

figure('position',[1091         235        1310        1026])

[mesh,data] = aplot.parse_mesh(mesh,i,data);
data.mesh   = mesh;

% convert network to aal90 space
%--------------------------------------------------------------------------
atlas = aplot.dotemplate('aal90');
rois  = aplot.get_roi_centres(atlas.template_sourcemodel.pos,atlas.all_roi_tissueindex);

atlas.template_sourcemodel.pos = rois;
atlas = rmfield(atlas,'all_roi_tissueindex');

% now register the real source model positions to the atlas ROIs
sm.pos = data.sourcemodel.pos;
reg    = aplot.interp_template(sm,rois);

% now apply the operator - M - to the connection matrix
NM     = reg.M;
NM     = NM/max(NM(:));          % rescale so not change amplitudes
[A,n]  = aplot.rw_edgenode(net); % read the matrix if necessary
S      = [min(A(:)) max(A(:))];  % get min/max
NL     = NM'*A*NM;               % apply M
A      = S(1) + ((S(2)-S(1))).*(NL - min(NL(:)))./(max(NL(:)) - min(NL(:)));
A(isnan(A)) = 0;

data.sourcemodel.pos = reg.pos;

% network  options
%--------------------------------------------------------------------------
netcmap = cmocean('balance'); % net color map
colbar  = 0;

data  = aplot.connections(data,A,colbar,i.write,i.fname,netcmap);