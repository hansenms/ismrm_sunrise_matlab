%%
%Clean up
close all;
clear all

noise_scale = 0.1;

%%
%Load Image & Sensitivity Maps
load im1.mat
load smaps_phantom.mat
load noise_covariances.mat

Rn = Rn_normal_8;
ncoils = size(smaps,3);
imsize = size(im1);
nx = imsize(1);
ny = imsize(2);
channel_im = smaps .* repmat(im1, [1 1 ncoils]);

data = ismrm_transform_image_to_kspace(channel_im, [1 2]);

x_lr = (1:32)+128-16;
y_lr = (1:32)+128-16;

f = hamming(32) * hamming(32)';
fmask = repmat(f, [1 1 ncoils]);
data_lr(x_lr, y_lr, :) = data(x_lr, y_lr, :) .* fmask;

im_lr = ismrm_transform_kspace_to_image(data_lr, [1 2]);

csm = ismrm_estimate_csm_mckenzie(im_lr);
ismrm_imshow(abs(csm));


%%
% Calibrate

csm = smaps;


unmix_sense = ismrm_calculate_sense_unmixing(acc_factor, csm, Rn);





%%
% 
x0 = 0;
y0 = 0;

data = zeros([imsize ncoils]);
for ic = 1:nc,
    data(:,:,ic) = smaps(x0, y0, ix) * ismrm_compute_gradient_encoding(-x0, -y0, nx, ny);
end
    
    
%%
% Reconstruct channel-by-channel images
channel_im_hat = ismrm_transform_kspace_to_image(data, [1 2]);

csm = smaps;

ccm = ismrm_compute_ccm(csm, Rn);




im_hat = sum( channel_im_hat .* conj(csm_walsh), 3 );

ismrm_imshow(abs(csm_walsh));
ismrm_imshow(angle(csm_walsh), []); colormap('hsv');

ismrm_imshow(abs(im_hat));
ismrm_imshow(angle(im_hat), [-pi pi]); colormap('hsv');

im_hat2 = sum( channel_im_hat, 3);


ismrm_imshow(abs(im_hat2));
ismrm_imshow(angle(im_hat2), [-pi pi]); colormap('hsv');
