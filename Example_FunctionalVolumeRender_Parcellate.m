
% - Render a functional (cortical) volume (nifti) on a template mesh
% - Compute a parcellation of the cortex (here, AAL-90)
% - Return mesh data, parcel volume & centroid (ROI) data / means.

% Load the parcel data: contains v (vertices) and vi (vector on integers
% describing which atlas/parcel/region each vertex belongs to
load LightAAL.mat

% The functional nifti volume to project
FunVol = 'Vis100B_Av.nii';

% The default mesh to render 
% (def1, 2 & 4 = cortical-only, 3=with cerebellum)
Surf = 'def2';

% Not specifying the 'method' defaults to ray-casting
afigure; D = atemplate('mesh',Surf,'overlay',FunVol,'post_parcel',{v vi});


% After rendering, all outputs are in struct 'D'.
%
% We can render the parcellation it returned:
% (I suggest using method 'spheres' for this)

afigure; atemplate('mesh',Surf,'overlay',D.post_parcel.data, ...
    'sourcemodel', D.post_parcel.pos , 'method', 'spheres' );


% We can also render the single atlas point means - i.e. for AAL90 we give
% it only 90 data points (D.post_parcel.ParcelMean). You can use the
% returned mean points for this:
%
% afigure; atemplate('mesh',Surf,'overlay',D.post_parcel.ParcelMean, ...
%    'sourcemodel', D.post_parcel.ParcelCent , 'method', 'spheres' );
%
%
% Since we have rendered AAL90 data, which there's already a default parcel
% mapping for, we can render the parcelation from the 90 points, using an
% additional method flag - this means we don't have to supploy the source
% locations / model:
%
afigure; atemplate('mesh',Surf,'overlay',D.post_parcel.data, ...
     'method', {'aal_light', 'spheres'} );
