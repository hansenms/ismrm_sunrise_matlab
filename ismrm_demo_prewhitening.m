%%
% Demo of effects of noise correlation

close all;
clear all;

%%
% Experiment Parameters
internal_cal_size = 20;
acc_factor = 2;
noise_scale = 0.1;

%%
%Load Data
fprintf('Loading Ground Truth Information\n');
load im1.mat
load smaps_phantom.mat
load noise_covariances.mat

Rn = Rn_broken_8;

csm = smaps;
im_shape = size(im1);
ncoils = size(smaps,3);
channel_im = smaps .* repmat(im1, [1 1 ncoils]);

%%
% create accelerated data
fprintf('Creating Accelerated Data\n');
noise_level = noise_scale*max(im1(:));
data = ismrm_transform_image_to_kspace(channel_im, [1 2]);

noise = ismrm_generate_correlated_noise(im_shape, Rn);
data = data + noise_level*noise;

sp = ismrm_generate_sampling_pattern(im_shape, acc_factor, internal_cal_size );
data_accel = data .* repmat((sp>0), [1 1 ncoils]);

im_full = ismrm_transform_kspace_to_image(data,[1,2]);
data_uniform_accel = data_accel .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]);
im_alias = ismrm_transform_kspace_to_image(data_uniform_accel,[1,2]);

num_reconstructions = 3;
titles = {'ignore correlation', 'prewhiten', 'use correlations'};
im_hat = zeros([im_shape num_reconstructions]);
gmap = zeros(size(im_hat));
noise_amplification = zeros(size(im_hat));

%%
% Reconstruction method 1: Ignore that we have correlated noise
fprintf('Recon, ignoring correlated noise\n');
unmix_sense1 = ismrm_calculate_sense_unmixing(acc_factor, csm) .* acc_factor;
im_hat(:,:,1) = sum(im_alias .* unmix_sense1, 3);

ccm = ismrm_compute_ccm(csm, Rn);
gmap(:,:,1) = ismrm_calculate_gmap(unmix_sense1, ccm, Rn, acc_factor);
noise_amplification(:,:,1) = ismrm_calculate_noise_amplification(unmix_sense1, Rn);
%%
% Reconstruction method 2: Apply a prewhitening step to create a set of
% virtual coils & data with no noise correlation and identical noise
% variance.
fprintf('Recon, apply prewhitening\n');
dmtx = ismrm_calculate_noise_decorrelation_mtx_from_covariance_mtx(Rn);
im_alias_prew = ismrm_apply_noise_decorrelation_mtx(im_alias, dmtx);
csm_prew = ismrm_apply_noise_decorrelation_mtx(csm, dmtx);
unmix_sense2 = ismrm_calculate_sense_unmixing(acc_factor, csm_prew) .* acc_factor;
im_hat(:,:,2) = sum(im_alias_prew .* unmix_sense2, 3);

ccm_prew = ismrm_compute_ccm(csm_prew);
gmap(:,:,2) = ismrm_calculate_gmap(unmix_sense2, ccm_prew, [], acc_factor);
noise_amplification(:,:,2) = sqrt(sum(abs(unmix_sense2).^2,3));


%%
% Reconstruction method 3: Pass the noise correlation matrix to the SENSE
% reconstruction.
fprintf('Recon, considering correlated noise\n');
unmix_sense3 = ismrm_calculate_sense_unmixing(acc_factor, csm, Rn) .* acc_factor;
im_hat(:,:,3) = sum(im_alias .* unmix_sense3, 3);

gmap(:,:,3) = ismrm_calculate_gmap(unmix_sense3, ccm, Rn, acc_factor);
noise_amplification(:,:,3) = ismrm_calculate_noise_amplification(unmix_sense3, Rn);

%%
% Plot results
im_max = max(abs(im1(:)));

ismrm_imshow(abs(im_hat), [0 1.0*im_max], [], titles);
ismrm_imshow(abs(im_hat-repmat(im1, [1 1 num_reconstructions])));
ismrm_imshow(gmap, [0 5]); colormap(jet); %colorbar;
ismrm_imshow(noise_amplification); colormap(jet); %colorbar;

