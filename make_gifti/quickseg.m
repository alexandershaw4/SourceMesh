function segmentedmri = quickseg(mri);

mri = ft_read_mri(mri);

template   = [fileparts(mfilename('fullpath')) '/T1.nii'];
fprintf('Warping MRI to template: %s\n',template);


cfg.coordsys = 'ctf';
mri = ft_volumerealign(cfg,mri,template);
cfg = [];

cfg.template   = template;
cfg.spmversion = 'spm8';
mri            = ft_volumenormalise(cfg,mri);
%mri            = ft_volumereslice([], mri);
cfg.brainsmooth= 6;
segmentedmri   = ft_volumesegment(cfg, mri);
