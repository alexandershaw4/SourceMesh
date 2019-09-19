% Example Script: 
%
% (1) Render an functional data (from volume) on a template brain, and
% (2) Project that 3D rendering onto a 2D plane

myvol = 'Vis100B_Av.nii';

afigure; D = atemplate('mesh','def2','overlay',myvol); % render the 3D mesh


% close, then get the mesh vertices, faces and colours for the 2D plot:
v  = D.mesh.vertices;
f  = D.mesh.faces;
c  = D.overlay.data;
cv = D.mesh.curvature;

afigure; [dv,uv,g] = flatten_mesh(v,f,cv,1); % make the awesome 2D plot