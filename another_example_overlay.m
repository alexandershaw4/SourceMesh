

figure,

load('example_aal90_data','M')
load DenseAAL.mat

D = atemplate('mesh','def5','sourcemodel',{v vi},'overlay',M);

% colormap(customcolormap_preset('pasteljet'))

D.mesh.h.LineStyle = 'none';
D.mesh.h.FaceAlpha=.4;
D.mesh.h.AmbientStrength=.8;

D.overlay.cb.CLim = [-2.5 3.5];