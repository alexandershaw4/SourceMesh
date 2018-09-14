
% Example fieldtrip segmentation usage with atemplate

load coregseg_0001.mat      % saved fieldtrip segmentation
 
segvol = segmentedmri.gray; % the gray matter volume


% extract surface & write gifti file
afigure,atemplate('mesh',segvol,'write','Alex')

% Now load up the saved gifti object & plot curvature overlay heatmap on
% leftside (for example)

afigure,atemplate('mesh','Alex.gii','overlay','curvature','hemi','l');
