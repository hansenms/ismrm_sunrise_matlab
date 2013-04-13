%%
% Demo of effects of noise correlation

close all;
clear all;

%%
% Experiment Parameters
internal_cal_size = 20;
acc_factor = 4;
noise_scale = 0.01;

%%
%Load Data
load im1.mat
%im1 = im1.';
im1 = phantom(256);
load smaps_phantom.mat
load noise_covariances.mat

Rn = Rn_normal_8;

csm = smaps;
im_shape = size(im1);
ncoils = size(smaps,3);
channel_im = smaps .* repmat(im1, [1 1 ncoils]);

%%
% create accelerated data
noise_level = noise_scale*max(im1(:));
data = ismrm_transform_image_to_kspace(channel_im, [1 2]);

noise = ismrm_generate_correlated_noise(im_shape, Rn);
data = data + noise_level*noise;

sp = ismrm_generate_sampling_pattern(im_shape, acc_factor, internal_cal_size );
data_accel = data .* repmat((sp>0), [1 1 ncoils]) .* acc_factor;

im_full = ismrm_transform_kspace_to_image(data,[1,2]);
data_uniform_accel = data_accel .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]);
im_alias = ismrm_transform_kspace_to_image(data_uniform_accel,[1,2]);% ./ acc_factor;

num_reconstructions = 3;
titles = {'mckenzie1', 'mckenzie2', 'walsh'};
im_hat = zeros([im_shape num_reconstructions]);
im_hat_full = zeros(size(im_hat));
gmap = zeros(size(im_hat));
noise_amplification = zeros(size(im_hat));

%%
% pull out cal data
cal_data = ismrm_extract_cal_data(data, sp);
cal_im = ismrm_transform_kspace_to_image(cal_data, [1,2], [im_shape ncoils]);
csm_mckenzie1 = ismrm_estimate_csm_mckenzie(cal_im);

f = hamming(size(cal_data,1)) * hamming(size(cal_data,2))';
fmask = repmat(f, [1 1 ncoils]);
im_lr = ismrm_transform_kspace_to_image(cal_data .* fmask, [1 2], [im_shape ncoils]);

csm_mckenzie2 = ismrm_estimate_csm_mckenzie(im_lr);
csm_walsh = ismrm_estimate_csm_walsh(im_lr);
 


%%
% Reconstruction method 1: mckenzie 1 maps
unmix_sense1 = ismrm_calculate_sense_unmixing(acc_factor, csm_mckenzie1, Rn);
im_hat(:,:,1) = sum(im_alias .* unmix_sense1, 3);

ccm = ismrm_compute_ccm(csm_mckenzie1, Rn);

im_hat_full(:,:,1) = sum(im_full .* ccm, 3);
gmap(:,:,1) = ismrm_calculate_gmap(unmix_sense1, ccm, Rn);
noise_amplification(:,:,1) = ismrm_calculate_noise_amplification(unmix_sense1, Rn);
%%
% mckenzie 2 maps
unmix_sense2 = ismrm_calculate_sense_unmixing(acc_factor, csm_mckenzie2, Rn, .1);
im_hat(:,:,2) = sum(im_alias .* unmix_sense2, 3);

ccm = ismrm_compute_ccm(csm_mckenzie2, Rn);

im_hat_full(:,:,2) = sum(im_full .* ccm, 3);
gmap(:,:,2) = ismrm_calculate_gmap(unmix_sense2, ccm, Rn);
noise_amplification(:,:,2) = ismrm_calculate_noise_amplification(unmix_sense2, Rn);

%%
% walsh maps
unmix_sense3 = ismrm_calculate_sense_unmixing(acc_factor, csm_walsh, Rn, .1);
im_hat(:,:,3) = sum(im_alias .* unmix_sense3, 3);

ccm = ismrm_compute_ccm(csm_walsh, Rn);

im_hat_full(:,:,3) = sum(im_full .* ccm, 3);
gmap(:,:,3) = ismrm_calculate_gmap(unmix_sense3, ccm, Rn);
noise_amplification(:,:,3) = ismrm_calculate_noise_amplification(unmix_sense3, Rn);

%%
% Plot results
im_max = max(abs(im_hat(:)));

ismrm_imshow(abs(im_hat), [0 0.3*im_max], [], titles);
ismrm_imshow(abs(im_hat_full));
ismrm_imshow(abs(im_hat-im_hat_full));
ismrm_imshow(gmap, [0 5]); colormap(jet); %colorbar;
ismrm_imshow(noise_amplification); colormap(jet); %colorbar;
