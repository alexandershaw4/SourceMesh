
% - Render a functional (cortical) volume (nifti) on a template mesh thresholded at top 20%
% - Inflate the mesh
% - Display the curvature of the brain in greyscale under the functional
% - Compute a parcellation of the cortex (here, AAL-90)
% - Return mesh data, parcel volume & centroid (ROI) data / means.
% - Render Right hemisphere only

% Load the parcel data: contains v (vertices) and vi (vector on integers
% describing which atlas/parcel/region each vertex belongs to
load LightAAL.mat

% The functional nifti volume to project
FunVol = 'Vis100B_Av.nii';

% The default mesh to render 
% (def1, 2 & 4 = cortical-only, 3=with cerebellum)
Surf = 'def1';

% the threshold
thr = 0.4;

% overlaying both curvature and function
afigure; D = atemplate('hemi','l','mesh',Surf,'inflate','overlay',...
			{'curvature',FunVol},'method','spheres',...
            'thresh',thr,'post_parcel',{v vi});


        
% use the output to render a1-value per parcel version
afigure;
atemplate('mesh','def1','inflate','overlay',D.post_parcel.ParVal,...
    'sourcemodel',D.post_parcel.pos,'method','spheres')