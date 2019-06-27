
% Plot functional overlay on structural image
%==========================================================================

% use template cortical + cerebellum mesh
struct = 'def3';

% use atlas-reduced code: aal116
funct = zeros(116,1);
funct(91:108)=-8:9;   % activations in cerebellum only

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

figure('position',[1091         235        1310        1026])

[mesh,data] = aplot.parse_mesh(mesh,i,data);
data.mesh   = mesh;

% functional overlay options
%--------------------------------------------------------------------------
data.overlay = [];
[y,data]     = aplot.parse_overlay(funct,data); % get functional vals & sourcemodel vertices

data.overlay.method  = {'aal116','spheres'}; % 'raycast', 'euclidean' or 'spheres'
data.overlay.peaks   = 0;
    
colbar = 0;
data   = aplot.overlay(data,y,i.write,i.fname,colbar);


