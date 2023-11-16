

figure,

load('example_aal90_data','M')
load DenseAAL.mat
load pastelmap

D = atemplate('mesh','def6','sourcemodel',{v vi},'overlay',M,'nocolbar');

colormap(map)

D.mesh.h.LineStyle = 'none';
D.mesh.h.FaceAlpha=.4;
D.mesh.h.AmbientStrength=.8;

D.overlay.cb.CLim = [-2.5 3.5];