%%
%Clean up
close all;
clear all

%%
% Parameters
noise_scale = 0.05;
kernel_shape = [5 7];

fprintf('Loading Ground Truth Information\n');
%%
%Load Image & Sensitivity Maps
load im1.mat
load smaps_phantom.mat
load noise_covariances.mat

ncoils = size(smaps,3);
Rn = Rn_normal_8;
im_shape = size(im1);
channel_im = smaps .* repmat(im1, [1 1 ncoils]);


%%
% create accelerated data
fprintf('Creating rFOV Accelerated Data\n');
noise = noise_scale * max(im1(:)) * ismrm_generate_correlated_noise(im_shape, Rn);

data = ismrm_transform_image_to_kspace(channel_im, [1 2]) + noise;

% halve FOV
data = data(:,1:2:end,:);
im_shape = [size(data,1) size(data,2)];

acc_factor = 2;
sp = ismrm_generate_sampling_pattern(im_shape, acc_factor, 20 );
data_accel = data .* repmat((sp>0), [1 1 ncoils]);

im_full = ismrm_transform_kspace_to_image(data,[1,2]);
data_uniform_accel = data_accel .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]);
im_alias = ismrm_transform_kspace_to_image(data_uniform_accel,[1,2]);

%%
% pull out cal data
cal_data = ismrm_extract_cal_data(data, sp);


%%
% Estimate sensitivitiy maps for channel combination
fprintf('Calibration - Estimating Coil Sensitivities\n');
f = hamming(size(cal_data,1)) * hamming(size(cal_data,2))';
fmask = repmat(f, [1 1 ncoils]);
im_lr = ismrm_transform_kspace_to_image(cal_data .* fmask, [1 2], [im_shape ncoils]);
csm_walsh = ismrm_estimate_csm_walsh(im_lr);
csm_mckenzie = ismrm_estimate_csm_mckenzie(im_lr);
csm_true = smaps(:,64 + (1:128),:);

csm_walsh = ismrm_normalize_shading_to_sos(csm_walsh);
csm_mckenzie = ismrm_normalize_shading_to_sos(csm_mckenzie);
csm_true = ismrm_normalize_shading_to_sos(csm_true);

ccm_true = ismrm_compute_ccm(csm_true);
ccm_mckenzie = ismrm_compute_ccm(csm_mckenzie);
ccm_walsh = ismrm_compute_ccm(csm_walsh);

%%
% Create unmixing images
fprintf('Calibration - Generate Unmixing Images\n');

jer_lookup_dd = ismrm_compute_jer_data_driven(cal_data, kernel_shape);
cal_im = ismrm_transform_kspace_to_image(cal_data, [1,2], 2 * size(cal_data));
jer_lookup_md = ismrm_compute_jer_model_driven(cal_im, kernel_shape);

num_recons = 5;
titles = {'SENSE true csm', 'SENSE walsh csm', 'SENSE mckenzie csm', 'PARS', 'GRAPPA'};
unmix = zeros([im_shape ncoils num_recons]);
unmix(:,:,:,1) = ismrm_calculate_sense_unmixing(acc_factor, csm_true, Rn) .* acc_factor;
unmix(:,:,:,2) = ismrm_calculate_sense_unmixing(acc_factor, csm_walsh, Rn) .* acc_factor;
unmix(:,:,:,3) = ismrm_calculate_sense_unmixing(acc_factor, csm_mckenzie, Rn) .* acc_factor;
unmix(:,:,:,4) = ismrm_calculate_jer_unmixing(jer_lookup_md, acc_factor, ccm_mckenzie, 0.001, false);
unmix(:,:,:,5) = ismrm_calculate_jer_unmixing(jer_lookup_dd, acc_factor, ccm_mckenzie, 0.001, false);

%%
% Perform reconstructions
fprintf('Perform reconstructions\n');
im_full   = abs(sum(im_full .* ccm_mckenzie, 3));

im_hat = zeros([im_shape, num_recons]);
im_diff = zeros([im_shape, num_recons]);


for recon_index = 1:num_recons,
    im_hat(:,:,recon_index) = abs(sum(im_alias .* unmix(:,:,:,recon_index), 3));
    im_diff(:,:,recon_index) = abs(im_hat(:,:,recon_index) - im_full);
end
    


ismrm_imshow(im_hat, [0 0.5 * max(abs(im_hat(:)))], [], titles);
ismrm_imshow(im_diff, [], [], titles);

ismrm_imshow(im_diff(:,:,4:5), [], [], {'|PARS-Full|', '|GRAPPA-Full|'})
