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

ismrm_imshow(abs(im1));
ismrm_imshow(angle(im1), [-pi pi]); colormap('hsv');

ismrm_imshow(abs(smaps));
ismrm_imshow(angle(smaps)); colormap('hsv');

ismrm_imshow(abs(channel_im));
ismrm_imshow(angle(channel_im)); colormap('hsv');


%%
% Synthesize Data
noise_scale = 1.0 * max(im1(:));
noise = noise_scale .* ismrm_generate_correlated_noise(imsize, eye(ncoils));
data = ismrm_transform_image_to_kspace(channel_im, [1 2]) + noise;


%%
% Reconstruct channel-by-channel images


channel_im_hat = ismrm_transform_kspace_to_image(data, [1 2]);

%csm_walsh = ismrm_estimate_csm_walsh(channel_im_hat);

csm_mckenzie = ismrm_estimate_csm_mckenzie(channel_im_hat);

im_hat = sum( channel_im_hat .* conj(csm_mckenzie), 3 );

%ismrm_imshow(abs(csm_walsh));
%ismrm_imshow(angle(csm_walsh), []);

ismrm_imshow(abs(im_hat));
ismrm_imshow(angle(im_hat), [-pi pi]); colormap('hsv');


