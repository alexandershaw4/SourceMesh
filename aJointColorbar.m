function aJointColorbar(D)
% Generate a colorbar when... 
% ... you've used atemplate to generate a mesh overlay that has both
% grayscale curvature and coloured function overlay values. It does this by
% generating a long colormap that is grayscale for the first 256 parts,
% then colour for the latter 256.
%
% The problem is that it means the default colourbar values are gibberish,
% so this function corrects them.
%
% Input is the structure outputted by atemplate.
%
% AS

% Find the axes
axes(D.mesh.h.Parent);
colorbar;

% the actual colorbar values
cbv = D.overlay.colbar_values;

% get the colorbar of the parent axes of the mesh plot
cb = get(D.mesh.h.Parent,'colorbar');

% flip the colour part of the bar
set(cb, 'YDir', 'reverse' );
set(cb, 'ylim', [1 256]);

% now generate some tickmarks ONLY for the color portion of the bar
points   = round(linspace(1,256,11));
cb.Ticks = points;
vals     = flipud( cbv(256+points) );
vals     = vals - vals(6);
vals(1:5)= -flipud(vals(7:end)); % correct rounding error

cb.TickLabels = vals;