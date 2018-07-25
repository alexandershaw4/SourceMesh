

% Overlay Examples:
%--------------------------------------------------------------------------



% Example (1): Functional overlay from the AAL90 atlas
%-----------------------------------------------------------------
figure('position',[186         290        1927         688])

overlay = randi([-6 6],90,1); % the functional values of length(90)

atemplate('overlay',overlay); % generate the image using defaults

exportbrains('MyOverlay')     % export high resolution images from all angles

slice3();                     % open a new figure showing left, right and top views



% Example (2): Functional overlay from a specified sourcemodel
%-----------------------------------------------------------------
figure('position',[186         290        1927         688])

load New_AALROI_6mm.mat       % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

overlay = randi([0 6],5061,1); % the functional values of length positions

atemplate('overlay',overlay,'sourcemodel', positions); % render the image



% Example (3): Project functional overlay from functional volume
%-----------------------------------------------------------------
figure('position',[186         290        1927         688])

myvolume = 'vol.nii';          % the functional volume filename

atemplate('overlay',myvolume); % render the image

% or myvolume could be a preloaded image volume, of dim 3.



% Example (4): Select the overlay registration/rendering method
%-----------------------------------------------------------------
figure('position',[186         290        1927         688])

load New_AALROI_6mm.mat        % load template sourcemodel

positions = template_sourcemodel.pos; % get positions

overlay = randi([0 6],5061,1);  % the functional values of length positions

method  = 'raycast';            % 'raycast', 'euclidean' or 'spheres' 

atemplate('overlay',overlay,... % render the image
          'sourcemodel', positions,...
          'method', method); 



% Example (4): Select the overlay registration/rendering method
%-----------------------------------------------------------------









