% This script provides examples of lots of different ways of using the plot
% function atemplate. Examples are structured by type with examples for
% each type of plot numbered.
%
% Sections:
% [A] - examples of how to plot meshes (templates, volumes etc).
% [B] - examples of how to plot functional overlays
% [C] - examples of how to plot networks
% [D] - examples of how to plot nodes
% [E] - examples of how to plot labels
% [F] - examples of how to make videos
% [G] - exports
%
%
% note: any combination of different types of plots is possible. e.g. you
% can specify to plot a network and functional overlay on a mesh you
% extracted from a volume.
%
% For more help, see the atemplate help.
%


% [A]: Basic Mesh Examples:
%--------------------------------------------------------------------------


% Example (1): Just use the default smoothed ICBM152 (81k)
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

atemplate();



% Example (2): Specify a gifti file (or structure)
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

g = 'MyMesh.gii';

atemplate('mesh', g);



% Example (3): Specify a nifti file (or volume)
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

g = 'MyVol.nii';

atemplate('mesh', g);




% [B]: Overlay Examples:
%--------------------------------------------------------------------------



% Example (1): Functional overlay from the AAL90 atlas
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

overlay = randi([-6 6],90,1); % the functional values of length(90)

atemplate('overlay',overlay); % generate the image using defaults

exportbrains('MyOverlay')     % export high resolution images from all angles

slice3();                     % open a new figure showing left, right and top views



% Example (2): Functional overlay from a specified sourcemodel
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

load New_AALROI_6mm.mat       % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

overlay = randi([0 6],5061,1); % the functional values of length positions

atemplate('overlay',overlay,'sourcemodel', positions); % render the image



% Example (3): Project functional overlay from functional volume
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

myvolume = 'vol.nii';          % the functional volume filename

atemplate('overlay',myvolume); % render the image

% or myvolume could be a preloaded image volume, of dim 3.

% option: Specify an affine transformation matrix for the functonal volume
% atemplate('overlay', myvolume, 'funcaffine', matrix);



% Example (4): Select the overlay registration/rendering method
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

load New_AALROI_6mm.mat        % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

overlay = randi([0 6],5061,1);  % the functional values of length positions

method  = 'raycast';            % 'raycast', 'euclidean' or 'spheres' 

atemplate('overlay',overlay,... % render the image
          'sourcemodel', positions,...
          'method', method); 

      
% option: when using ray casting, specify a depth vector (positions along
% the face normal line) to search:
% atemplate( ... , 'method', 'raycast', 'depth', -1.5:0.05:1.5);


% Example (5): Transform an overlay into AAL atlas space (roughly)
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

load New_AALROI_6mm.mat        % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

overlay = randi([0 6],5061,1);  % the functional values of length positions

atlas   = 'aal90';              % atlas to use: choose aal90 / aal78 / aal58

atemplate('overlay',overlay,'sourcemodel',positions,'template','aal90');


% Example (5b): Transform an overlay into AAL atlas space (roughly) AND
%               then use an existing mapping from 90 regions to template
%               [MNI] surface parcelation:
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

load New_AALROI_6mm.mat        % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

overlay = randi([0 6],5061,1);  % the functional values of length positions

atlas   = 'aal90';              % atlas to use: choose aal90 / aal78 / aal58

atemplate('overlay',overlay,'sourcemodel',positions,'template','aal90',...
        'method',{'aal_light','spheres'});

    
% Example (5c): as above, but with a specific ROI specified
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

load New_AALROI_6mm.mat        % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

overlay = randi([0 6],5061,1);  % the functional values of length positions

atlas   = 'aal90';              % atlas to use: choose aal90 / aal78 / aal58

roi     = {'Frontal_Sup_Orb_L' 'Supp_Motor_Area_L' 'Heschl_R'};

atemplate('mesh','def4','overlay',overlay,'sourcemodel',positions,'template','aal90',...
        'method',{'aal_light','spheres'},'roi',roi);


% Example (6): Perform a spatial pca or find local maxima 
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

load New_AALROI_6mm.mat         % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

overlay = randi([0 6],5061,1);  % the functional values of length positions

atemplate('overlay',overlay,'sourcemodel',positions,'pca');        % for pca
atemplate('overlay',overlay,'sourcemodel',positions,'components'); % for lm
atemplate('overlay',overlay,'sourcemodel',positions,'peaks');  % for peak labels (when using AAL90)



% [C]: Network Examples:
%--------------------------------------------------------------------------


% Example (1): Plot a network from the AAL90 atlas on default mesh
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

net                                        = zeros(90,90);
net(randi([1 90],20,1),randi([1 90],20,1)) = randi([-12 16],20,20);

atemplate('network', net);



% Example (2): Plot a network from a specified source model
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

load New_AALROI_6mm.mat         % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

nv                       = length(positions);
net                      = zeros(nv,nv);
net([1 4 78],[5 6 7 19]) = randi([-6 6],3,4);

atemplate('network', net, 'sourcemodel', positions);

% Specify a colormap for this:
atemplate('network', net, 'sourcemodel', positions, 'netcmap', jet)



% Example (3): Plot a network from .node & .edge files:
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

netfile = 'Network.edge';      % specify the .edge (it finds the .node)

atemplate('network', netfile); % render the network



% Example (4): Transform network into AAL atlas space (roughly)
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

netfile = 'Network.edge';      % specify the .edge (it finds the .node)

atlas   = 'aal90';             % atlas to use: choose aal90 / aal78 / aal58

atemplate('overlay',netfile,'template','aal90');



% [D]: Node Examples:
%--------------------------------------------------------------------------


% Example (1): Plot selected nodes of the AAL atlas
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

node = randi([0 1],90,1); % random bunch of the nodes to display

atemplate('nodes',node)   % render them on the default mesh


% Example (2): Plot selected nodes of a specified source model
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

load New_AALROI_6mm.mat         % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

node = randi([0 3],length(positions),1); % random bunch of the nodes to display

atemplate('nodes',node)   % render them on the default mesh


% [E]: Labels Examples:
%--------------------------------------------------------------------------


% Example (1): Plot all of the AAL90 labels on the default mesh
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

atemplate('labels')


% Example (2): Plot the AAL90 labels of the plotted nodes only:
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

node = randi([0 1],90,1); % random bunch of the nodes to display

atemplate('nodes',node,'labels')


% Example (3): Plot the AAL90 labels of the plotted network nodes
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

net                                        = zeros(90,90);
net(randi([1 90],20,1),randi([1 90],20,1)) = randi([-12 16],20,20);

atemplate('network', net, 'labels');


% Example (4): Add labels to a custom source model at the centre of the
% voxels belonging to that ROI:
%-----------------------------------------------------------------
figure('position',[481         176        1204         987])

load New_AALROI_6mm.mat               % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

% a 1-by-num-vertices vector containing indices of the roi this vertex belongs to
all_roi_tissue = zeros(length(positions),1);
all_roi_tissue(1   :1000) = 1; % the first 300 points represent ROI 1
all_roi_tissue(2000:3000) = 2; % the next 300 points represent  ROI 2
all_roi_tissue(4000:5000) = 3; % the next 300 points represent  ROI 3

% the labels for each roi
labels         = {'ROI 1' 'ROI 2' 'ROI 3'}; 

atemplate('sourcemodel',positions,'labels', all_roi_tissue, labels); 





