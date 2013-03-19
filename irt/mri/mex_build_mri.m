% mex_build_mri
% run matlab's "mex" command to "compile" the mri-related code
% into mex files.
% only users on unsupported systems, e.g., PCs, will need to do this

dir_current = pwd;
dir_mri = path_find_dir('mri');
cd(dir_mri)
mex exp_xform_mex.c
mex mri_exp_mult_mex.c
cd(dir_current)
