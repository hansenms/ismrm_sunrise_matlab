%%
%Clean up
close all;
clear all

%%
%Load Image & Sensitivity Maps
load ../im1.mat
load ../smaps_phantom.mat
im1 = im1;
ncoils = size(smaps,3);
imsize = size(im1);

channel_im = smaps .* repmat(im1, [1 1 ncoils]);



%%
% create accelerated data
close all;
acc_factor = 4;
noise_level = 0.05*max(im1(:));
[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 20);
data_noise = data + noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);

data_noise = data_noise .* acc_factor;

img_alias_noise = ismrm_transform_kspace_to_image(data_noise .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]),[1,2]);

%%
% pull out cal data
ky_projection = max(sp, [], 1);
kx_projection = max(sp, [], 2);
cal_indx = find(kx_projection > 1);
cal_indy = find(ky_projection > 1);

kx_cal_bounds = [min(cal_indx):max(cal_indx)];
ky_cal_bounds = [min(cal_indy):max(cal_indy)];
cal_data = data_noise(kx_cal_bounds, ky_cal_bounds, :);
nx = size(im1,1);
ny = size(im1,2);

%%
% 
csm = smaps;
ccm = ismrm_compute_ccm(csm) .* (1.0/acc_factor);

%ccm = zeros(size(csm));
%ccm(:,:,1) = 1;

kernel_shape = [5 15];

fprintf('computing jer lookup\n')
tic
%jer_lookup = compute_jer_model_driven(channel_im, kernel_shape);

jer_lookup = compute_jer_data_driven(cal_data, kernel_shape);

%cal_im = ismrm_transform_kspace_to_image(cal_data, [1,2], 2*size(cal_data));
%jer_lookup = compute_jer_model_driven(cal_im, kernel_shape);

toc

unmix = ismrm_calculate_jer_unmixing(jer_lookup, acc_factor, ccm, true);
ismrm_imshow(abs(unmix));
ismrm_imshow(abs(sum(img_alias_noise .* unmix,3))); 
