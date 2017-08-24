

%IN = Vol2SurfAS('mymri.mri','ctf',.15);

F = figure;
h1 = subplot(121);
atemplate('gifti',IN,'fighnd',h1,'network',edges,'labels','nocolbar')
view(0,0);

h2 = subplot(122);
atemplate('gifti',IN,'fighnd',h2,'overlay',t_ol,'nocolbar');
view(0,0);

set(gcf, 'Position', [100, 100, 2400, 1000])

im2vid_multi(F,'test2');