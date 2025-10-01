

afigure;

load('example_aal90_data','M')
load DenseAAL.mat
load pastelmap

%D = atemplate('mesh','def6','sourcemodel',{v vi},'overlay',M,'nocolbar');


opts = struct('normal_depth_mm',2.5,'normal_steps',7,'aggregate','max', ...
              'inpaint_iters',1,'smooth_iters',3);

D = atemplate('mesh','def6', ...
              'sourcemodel',{v vi},...
              'overlay',M, ...
              'overlay_opts',opts, ...
              'method','raycast');
axis equal;

colormap(map)

D.mesh.h.LineStyle = 'none';
D.mesh.h.FaceAlpha=.4;
D.mesh.h.AmbientStrength=.8;

D.overlay.cb.CLim = [-2.5 3.5];