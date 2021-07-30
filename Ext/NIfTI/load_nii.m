%  Load NIFTI dataset. Support both *.nii and *.hdr/*.img file extension.
%  If file extension is not provided, *.hdr/*.img will be used as default.
%
%  A subset of NIFTI sform/qform transform is performed in load_nii
%  to translate, flip, rotate (limited to N*90 degree), intensity 
%  scaling etc. Other transforms (any degree rotation, shears, etc.)
%  are not supported, because in those transforms, each voxel has 
%  to be repositioned, interpolated, and whole image(s) will have 
%  to be reconstructed. If an input data (nii) can not be handled,
%  The program will exit with an error message "Transform of this 
%  NIFTI data is not supported by the program". After the transform,
%  nii will be in RAS orientation, i.e. X axis from Left to Right,
%  Y axis from Posterior to Anterior, and Z axis from Inferior to
%  Superior. The RAS orientation system sometimes is also referred
%  as right-hand coordinate system, or Neurologist preferred system.
%  
%  Usage: [nii] = load_nii(filename, [img_idx], [old_RGB], [preferredForm])
%  
%  filename - NIFTI file name.
%  
%  img_idx    - a numerical array of image indices. Only the specified
%	images will be loaded. If there is no img_idx, all available 
%	images will be loaded.
%
%  old_RGB    - a boolean variable to tell difference from new RGB24 from old
%       RGB24. New RGB24 uses RGB triple sequentially for each voxel, like
%       [R1 G1 B1 R2 G2 B2 ...]. Analyze 6.0 developed by AnalyzeDirect uses
%       old RGB24, in a way like [R1 R2 ... G1 G2 ... B1 B2 ...] for each
%       slices. If the image that you view is garbled, try to set old_RGB
%       variable to 1 and try again, because it could be in old RGB24.
%
%  preferredForm  - selects which transformation from voxels to RAS
%	coordinates; values are s,q,S,Q.  Lower case s,q indicate
%	"prefer sform or qform, but use others if preferred not present". 
%	Upper case indicate the program is forced to use the specificied
%	tranform or fail loading.  Default is 's'.	- Jeff Gunter
%
%  The number of images scans can be obtained from get_nii_frame, or
%  simply: hdr.dime.dim(5)
%  
%  Returned values:
%  
%  nii.hdr - struct with NIFTI header fields.
%  nii.filetype - Analyze format (0); NIFTI .hdr/.img (1); NIFTI .nii (2)
%  nii.fileprefix - NIFTI filename without extension.
%  nii.machine - machine string variable.
%  nii.img_idx - Indices of images to be loaded.
%  nii.img - 3D (or 4D) matrix of NIFTI data.
%  
%  Part of this file is copied and modified under GNU license from
%  MRI_TOOLBOX developed by CNSP in Flinders University, Australia
%  
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%  
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
function [nii] = load_nii(filename, img_idx, old_RGB, preferredForm)

   if ~exist('filename','var'),
      error('Usage: [nii] = load_nii(filename, [img_idx], [old_RGB])');
   end
   
   if ~exist('img_idx','var'), img_idx = []; end
   if ~exist('old_RGB','var'), old_RGB = 0; end
   if ~exist('preferredForm','var'), preferredForm= 's'; end     % Jeff
   
   %  Read the dataset header
   %
   [nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr(filename);

   %  Read the header extension
   %
%   nii.ext = load_nii_ext(filename);
   
   %  Read the dataset body
   %
   [nii.img,nii.hdr] = ...
	load_nii_img(nii.hdr,nii.filetype,nii.fileprefix,nii.machine,img_idx,old_RGB);
   
   %  Perform some of sform/qform transform
   %
      
   nii = xform_nii(nii, preferredForm);

   return					% load_nii

