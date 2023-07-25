% script to make plots using patch when using a triangulated source model
% and functional overlay values

p = patch('faces',HCP360.f,'vertices',HCP360.v);
p.FaceVertexCData = o;
shading interp;
camlight left; camlight right;
lighting gouraud;
material dull
view(3)