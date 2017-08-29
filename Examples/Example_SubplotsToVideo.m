

G  = read_nv;
iG = spm_mesh_inflate(G);
iG = gifti(iG);

F  = figure;
h1 = subplot(121);
atemplate('gifti',iG,'fighnd',h1,'network',randi([0 4],90,1),'labels','nocolbar')
view(0,0);

h2 = subplot(122);
atemplate('gifti',iG,'fighnd',h2,'overlay',randi([0 4],90,1),'nocolbar');
view(0,0);

set(gcf, 'Position', [100, 100, 2400, 1000])

im2vid_multi(F,'test3');