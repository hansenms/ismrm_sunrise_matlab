%%
%Clean up
close all;
clear all

acc_factor = 2;
cal_noise_scale = 0.05;
accel_noise_scale = 0.05;
cal_shape = [64 6];
kernel_shape = [5 3];

%%
%Load Image & Sensitivity Maps
fprintf('Loading Ground Truth Information\n');
load im1.mat
load smaps_phantom.mat
load noise_covariances.mat
%smaps = smaps(:,:,1:2:end);
im1 = im1.';

ncoils = size(smaps,3);
Rn = eye(ncoils);
imsize = size(im1);
nx = imsize(1);
ny = imsize(2);

pixel_mask = sum(abs(smaps),3) > 0;

%%
% Create Calibration Data
fprintf('Calibrating - Creating Calibration Data\n');

channel_im = smaps .* repmat(im1, [1 1 ncoils]);
cal_data = ismrm_transform_image_to_kspace(channel_im, [1 2], cal_shape);
noise = cal_noise_scale * max(im1(:)) * ismrm_generate_correlated_noise(cal_shape, Rn);
cal_data = cal_data + noise;

f = hamming(cal_shape(1)) * hamming(cal_shape(2))';
fmask = repmat(f, [1 1 ncoils]);
filtered_cal_data = cal_data .* fmask;

cal_im = ismrm_transform_kspace_to_image(filtered_cal_data, [1 2], imsize);
cal_im = cal_im .* repmat(pixel_mask, [1 1 ncoils]);

fprintf('Calibrating - Estimating Coil Sensitivities\n');
csm_walsh = ismrm_estimate_csm_walsh(cal_im);
csm_mckenzie = ismrm_estimate_csm_mckenzie(cal_im);

csm_walsh = ismrm_normalize_shading_to_sos(csm_walsh);
csm_mckenzie = ismrm_normalize_shading_to_sos(csm_mckenzie);
csm_true = ismrm_normalize_shading_to_sos(smaps);

ccm_true = ismrm_compute_ccm(csm_true);
ccm_mckenzie = ismrm_compute_ccm(csm_mckenzie);
ccm_walsh = ismrm_compute_ccm(csm_walsh);

%%
% Calibrate
fprintf('Calibrating - Generating Unmixing Images\n');
jer_lookup_dd = ismrm_compute_jer_data_driven(cal_data, kernel_shape);
jer_lookup_md1 = ismrm_compute_jer_model_driven(cal_im, kernel_shape);

unfiltered_cal_im = ismrm_transform_kspace_to_image(cal_data, [1,2], 2 * size(cal_data));
jer_lookup_md2 = ismrm_compute_jer_model_driven(unfiltered_cal_im, kernel_shape);

num_recons = 6;
titles = {'SENSE true csm', 'SENSE mckenzie csm', 'SENSE walsh csm', 'PARS filtered cal im', 'PARS unfiltered cal im', 'GRAPPA'};

unmix = zeros([imsize ncoils num_recons]);
unmix(:,:,:,1) = ismrm_calculate_sense_unmixing(acc_factor, csm_true, Rn, 0) .* acc_factor;
unmix(:,:,:,2) = ismrm_calculate_sense_unmixing(acc_factor, csm_mckenzie, Rn, 0) .* acc_factor;
unmix(:,:,:,3) = ismrm_calculate_sense_unmixing(acc_factor, csm_walsh, Rn, 0) .* acc_factor;
unmix(:,:,:,4) = ismrm_calculate_jer_unmixing(jer_lookup_md1, acc_factor, ccm_mckenzie, 0, false);
unmix(:,:,:,5) = ismrm_calculate_jer_unmixing(jer_lookup_md2, acc_factor, ccm_mckenzie, 0, false);
unmix(:,:,:,6) = ismrm_calculate_jer_unmixing(jer_lookup_dd, acc_factor, ccm_mckenzie, 0, false);


%%
% Create Accelerated Data
fprintf('Creating Accelerated Data\n');
noise = accel_noise_scale * max(im1(:)) * ismrm_generate_correlated_noise(imsize, Rn);
data = ismrm_transform_image_to_kspace(channel_im, [1 2]) + noise;
sp = ismrm_generate_sampling_pattern(size(im1), acc_factor, 0);
data_accel = data .* repmat(sp == 1 | sp == 3,[1 1 ncoils]);
im_alias = ismrm_transform_kspace_to_image(data_accel,[1,2]);% ./ acc_factor;

im_full = ismrm_transform_kspace_to_image(data, [1, 2]);

%ccm_true = ismrm_compute_ccm(smaps, Rn);
im_full = abs(sum(im_full .* ccm_true, 3));

%%
% Analyze Reconstruction Candidates
fprintf('Performing and Analyzing Reconstruction Candidates\n');
aem = zeros([imsize num_recons]);
gmap = zeros([imsize, num_recons]);
im_hat = zeros([imsize, num_recons]);
im_diff = zeros([imsize, num_recons]);

signal_mask = imclose(im1>100.0, strel('disk', 5)); ismrm_imshow(signal_mask, [0 1]);


for recon_index = 1:num_recons,
    aem(:,:,recon_index) = ismrm_calculate_aem(signal_mask, csm_true, unmix(:,:,:,recon_index), acc_factor);
    gmap(:,:,recon_index) = ismrm_calculate_gmap(unmix(:,:,:,recon_index), ccm_true, Rn, acc_factor);
    im_hat(:,:,recon_index) = abs(sum(im_alias .* unmix(:,:,:,recon_index), 3));
    im_diff(:,:,recon_index) = abs(im_hat(:,:,recon_index) - im_full);
end
    


ismrm_imshow(aem, [0 0.1], [], titles); colormap(jet);
ismrm_imshow(gmap, [0 6], [], titles); colormap(jet);
ismrm_imshow(im_hat, [0 0.7 .* max(im_hat(:))], [], titles);
ismrm_imshow(im_diff, [0 0.1 * max(im_hat(:))], [], titles);