function linksubplots(ha)
% Link subplot axes for roation purposes.
% ha = figure handle 
%
% If fig reopened from .fig do:
%   h  = figopen('blah.fig')
%   ha = allchild(h); 
%   linksubplots(ha)
%
% AS

Link = linkprop(ha, ...
       {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);