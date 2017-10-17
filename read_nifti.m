function D = read_nifti(filein)

N = spm_vol(filein);
D = spm_read_vols(N);