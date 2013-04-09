%%
%Clean up
close all;
clear all

%%
%Load Image & Sensitivity Maps
load ../im1.mat
load ../smaps_phantom.mat
ncoils = size(smaps,3);
imsize = size(im1);

channel_im = smaps .* repmat(im1, [1 1 ncoils]);

%ismrm_imshow(abs(im1));
%ismrm_imshow(angle(im1), [-pi pi]); colormap('hsv');

%ismrm_imshow(abs(smaps));
%ismrm_imshow(angle(smaps)); colormap('hsv');

%ismrm_imshow(abs(channel_im));
%ismrm_imshow(angle(channel_im)); colormap('hsv');

%ismrm_imshow(abs(sum(channel_im, 3)));

%%
%
close all;
acc_factor = 4;
noise_level = 0.05*max(im1(:));
[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 64);
data_noise = data + noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);

%data       = data * acc_factor;
data_noise = data_noise .* acc_factor;

%img_alias = ismrm_transform_kspace_to_image(data .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]),[1,2]);
img_alias_noise = ismrm_transform_kspace_to_image(data_noise .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]),[1,2]);
%ismrm_imshow(abs(img_alias));
%ismrm_imshow(abs(img_alias_noise));

% pull out cal data
kx_projection = max(sp, [], 1);
ky_projection = max(sp, [], 2);
cal_indx = find(kx_projection > 1);
cal_indy = find(ky_projection > 1);

kx_cal_bounds = [min(cal_indx):max(cal_indx)];
ky_cal_bounds = [min(cal_indy):max(cal_indy)];
cal_data = data_noise(kx_cal_bounds, ky_cal_bounds, :);
cal_data = data(kx_cal_bounds, ky_cal_bounds, :);

nx = size(im1,1);
ny = size(im1,2);

%ismrm_imshow(abs(channel_im));

cal_data2 = ismrm_transform_image_to_kspace(channel_im, [1 2]);
unmix_grappa = ismrm_calculate_grappa_unmixing2(cal_data2, [3 2], acc_factor, [], [], [], [nx ny], true);


%ismrm_imshow(abs(sum(img_alias .* unmix_grappa,3))); colorbar; title('GRAPPA (noiseless)');

ismrm_imshow(abs(sum(img_alias_noise .* unmix_grappa,3))); %title('GRAPPA (noise)');
%ismrm_imshow(gmap_grappa); colorbar; title('GRAPPA g-factor');

%showimage(sum(img_alias_noise .* unmix_grappa,3),[2 2 3]); colorbar; title('GRAPPA');
%showimage(gmap_grappa,[2 2 4]); colorbar; title('GRAPPA g-factor');
%colormap(gray);
%set(gcf,'color','w');