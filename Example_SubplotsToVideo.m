

F = figure;
h1 = subplot(121);
atemplate('gifti',IN,'fighnd',h1,'labes')
view(0,0);

h2 = subplot(122);
atemplate('gifti',IN,'fighnd',h2,'overlay',randi([0 4],90,1));
view(0,0);

set(gcf, 'Position', [100, 100, 2400, 1000])

im2vid_multi(F,'test2');