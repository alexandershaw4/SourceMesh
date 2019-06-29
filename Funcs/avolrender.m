function D = avolrender(filename)
% a quick-to-execute wrapper on atemplate.m
%
%
%

% assume we want to register to the cortical AAL parcellation
load DenseAAL.mat

% The functional nifti volume to project
FunVol = filename;

% The default mesh to render 
% (def1, 2 & 4 = cortical-only, 3=with cerebellum)
Surf = 'def2';

% Not specifying the 'method' defaults to ray-casting
afigure; D = atemplate('mesh',Surf,'overlay',FunVol,'post_parcel',{v vi});
