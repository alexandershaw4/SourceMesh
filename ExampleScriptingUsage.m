
% Plot functional overlay on structural image
%==========================================================================

% from nifti's:
struct = 'GreyMRI.nii';
funct  = 'nTEST_Gamma,60-90,60-90Hz_pairedt.nii';

% or from gifti's:
struct = 'TestMesh.gii';
funct  = 'TestMeshOverlay.gii';


% make struct surface
%--------------------------------------------------------------------------
i.g  = struct;
data = aplot.sort_sourcemodel([],i);
[mesh,data] = aplot.get_mesh(i,data);       

% structural options
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

[mesh,data] = aplot.parse_mesh(mesh,i,data);
data.mesh   = mesh;

% functional overlay options
%--------------------------------------------------------------------------
data.overlay = [];
[y,data]     = aplot.parse_overlay(funct,data); % get functional vals & sourcemodel vertices

data.overlay.method  = 'raycast'; % 'raycast', 'euclidean' or 'spheres'
data.overlay.peaks   = 0;
    
colbar = 0;
data   = aplot.overlay(data,y,i.write,i.fname,colbar);


