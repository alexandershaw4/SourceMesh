function CompareProjectionMethods(sourcemodel,overlay)

% Compare the 3 overlay methods

figure('position',[186         290        1927         688])

% Ray casting
s(1) = subplot(131);
atemplate('overlay',overlay,...
          'sourcemodel',sourcemodel,...
          'method','raycast',...
          'fighnd', s(1))
title('Ray casting');

% Euclidean search
s(2) = subplot(132);
atemplate('overlay',overlay,...
          'sourcemodel',sourcemodel,...
          'method','euclidean',...
          'fighnd', s(2))
title('Euclidean search');

% Inflated sphere / trap radius
s(3) = subplot(133);
atemplate('overlay',overlay,...
          'sourcemodel',sourcemodel,...
          'method','spheres',...
          'fighnd', s(3))
title('Inflated spherical trap radius');


linksubplots(s);

