%%
%Load Data
close all;
clear all;

load im1.mat
load smaps_phantom.mat
load noise_covariances.mat

L_normal_8 = chol(Rn_normal_8,'lower');
L_broken_8 = chol(Rn_broken_8,'lower');
ncoils = size(smaps,3);

%Some settings
acc_factor = 4;
noise_level = 0.05*max(im1(:));

%Sample Cartesian Data
noise_white = noise_level*complex(randn(size(im1,1),size(im1,2),size(smaps,3)),randn(size(im1,1),size(im1,2),size(smaps,3)));
noise_color = reshape(permute(L_broken_8 * permute(reshape(noise_white, numel(noise_white)/ncoils,ncoils),[2 1]),[2 1]),size(noise_white));


[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 0);
data = data + noise_color .* repmat(sp > 0,[1 1 ncoils]);


%Save the variables that we need
save hansen_exercises.mat noise_color data sp smaps