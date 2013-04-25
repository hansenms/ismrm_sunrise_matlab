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

%Non-Cartesian data
[k_spiral,w_spiral] = ismrm_calculate_spiral_trajectory(acc_factor);

%Prepare NUFFT
N = [size(im1,1) size(im1,2)];
J = [5 5];
K = N*2;
nufft_st = nufft_init(k_spiral*2*pi,N,J,K,N/2,'minmax:kb');

%Sample
data_noiseless = nufft(repmat(im1,[1 1 size(smaps,3)]).*smaps,nufft_st)  ./ sqrt(prod(N));
noise_spiral = noise_level*complex(randn(size(data_noiseless)),randn(size(data_noiseless)));
data_spiral = data_noiseless + noise_spiral;

%Finally create a blurred regularization image:
reg_size = 32;
reg_img = ismrm_transform_kspace_to_image(padarray(hamming(reg_size)*hamming(reg_size)',[(size(im1,1)-reg_size)/2 (size(im1,2)-reg_size)/2],0,'both').*ismrm_transform_image_to_kspace(im1));
reg_img = abs(reg_img)+1;

%Save the variables that we need
save -v6 hansen_exercises.mat noise_color data sp smaps data_spiral noise_spiral k_spiral w_spiral reg_img