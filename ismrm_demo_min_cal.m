%%
%Clean up
close all;
clear all

num_cal_lines = 32;
acc_factor = 4;
kernel_shape = [5 7];


%%
%Load Image & Sensitivity Maps
load ../im1.mat
load ../smaps_phantom.mat
%smaps = smaps(:,:,1:2:end);
ncoils = size(smaps,3);
channel_im = smaps .* repmat(im1, [1 1 ncoils]);


%%
% create accelerated data
noise_level = 0.01*max(im1(:));
data = ismrm_transform_image_to_kspace(channel_im, [1 2]);
data = data + noise_level*complex(randn(size(data)),randn(size(data)));
% halve FOV
%data = data(:,1:2:end,:);
im_shape = [size(data,1) size(data,2)];


sp = ismrm_generate_sampling_pattern(im_shape, acc_factor, num_cal_lines);
data_accel = data .* repmat((sp>0), [1 1 ncoils]) .* acc_factor;

im_full = ismrm_transform_kspace_to_image(data,[1,2]);
data_uniform_accel = data_accel .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]);
im_alias = ismrm_transform_kspace_to_image(data_uniform_accel,[1,2]) ./ acc_factor;

%%
% pull out cal data
cal_data = ismrm_extract_cal_data(data, sp);


%%
% Estimate sensitivitiy maps for channel combination
f = hamming(size(cal_data,1)) * hamming(size(cal_data,2))';
fmask = repmat(f, [1 1 ncoils]);
im_lr = ismrm_transform_kspace_to_image(cal_data .* fmask, [1 2], [im_shape ncoils]);
csm = ismrm_estimate_csm_walsh(im_lr);

csm = smaps;
ccm = ismrm_compute_ccm(csm) .* (1.0/acc_factor);

%%

fprintf('computing jer lookup\n')


jer_lookup_dd = compute_jer_data_driven(cal_data, kernel_shape);
%jer_lookup_dd = jer_ref(cal_data, kernel_shape);
cal_im = ismrm_transform_kspace_to_image(cal_data, [1,2], 2 * size(cal_data));
jer_lookup_md = compute_jer_model_driven(cal_im, kernel_shape);



unmix_grappa = ismrm_calculate_jer_unmixing(jer_lookup_dd, acc_factor, ccm, true);
unmix_pars   = ismrm_calculate_jer_unmixing(jer_lookup_md, acc_factor, ccm, true);
unmix_sense = ismrm_calculate_sense_unmixing(acc_factor, csm);
%unmix_grappa2 = ismrm_calculate_grappa_unmixing(data_accel, [5 4], acc_factor, (sp>1), csm);

num_images = 4;
m_array = zeros([im_shape num_images]);
m_array(:,:,1)   = abs(sum(im_full .* ccm, 3));
m_array(:,:,2) = abs(sum(im_alias .* unmix_grappa,3));
m_array(:,:,3) = abs(sum(im_alias .* unmix_pars,3));
m_array(:,:,4) = abs(sum(im_alias .* unmix_sense,3));
%m_array(:,:,5) = abs(sum(im_alias .* unmix_grappa2,3));

titles = {'unaccelerated', 'grappa', 'pars', 'sense'};
ismrm_imshow(m_array, [], [2 2], titles);

diff_array = abs(m_array - repmat(m_array(:,:,1), [1 1 num_images]));
ismrm_imshow(diff_array, [], [], titles);
%ismrm_imshow(m_grappa);
%ismrm_imshow(m_pars);
%ismrm_imshow(m_sense);
%ismrm_imshow(m_full);

%ismrm_imshow(abs(m_full - m_grappa));
%ismrm_imshow(abs(m_full - m_pars));
%ismrm_imshow(abs(m_full - m_sense));


