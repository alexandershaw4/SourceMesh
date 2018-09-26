function CompareProjectionMethods(sourcemodel,overlay)

% Compare the 3 overlay methods

figure('position',[186         290        1927         688])

mesh = 'def1';

% Ray casting
s(1) = subplot(131);
atemplate('mesh',mesh,'overlay',overlay,...
          'sourcemodel',sourcemodel,...
          'method','raycast',...
          'fighnd', s(1),'hemi','l')
title('Ray casting');

% Euclidean search
s(2) = subplot(132);
atemplate('mesh',mesh,'overlay',overlay,...
          'sourcemodel',sourcemodel,...
          'method','euclidean',...
          'fighnd', s(2),'hemi','l')
title('Euclidean search');

% Inflated sphere / trap radius
s(3) = subplot(133);
atemplate('mesh',mesh,'overlay',overlay,...
          'sourcemodel',sourcemodel,...
          'method','spheres',...
          'fighnd', s(3),'hemi','l')
title('Inflated spherical trap radius');


linksubplots(s);

