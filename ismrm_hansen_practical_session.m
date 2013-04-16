%% ISMRM Sunrise Practical Session 
%
% This document contains the second set of practical exercises for the
% ISMRM course on parallel imaging.

%% Excercise Data
% All the data used in this set of exercises can be found in the file
% hansen_exercises.mat.
% We will start by clearing the workspace and loading the data

close all; clear all;
load hansen_exercises.mat
whos

%% Noise Pre-Whitening
%
%  The purpose of this exercise is to see the effects of noise
%  pre-whitening. We will use a SENSE reconstruction as an example
%
% $$ \tilde{\mathbf{x}} =  \arg \min_{\mathbf{x}} \left\{ \left\| \mathbf{E}\mathbf{x}-\mathbf{m} \right\|_2 \right\} $$
%
% where $\mathbf{x}$ is the reconstructed image, $\mathbf{E}$ is the
% encoding (sensitivity matrix) and $\mathbf{s}$ is the measured data.
%
% Let's start with a naive SENSE reconstruction

smask = (sp == 1 | sp == 3);
ncoils = size(smaps,3);
nelements = numel(smaps)/ncoils;

acc_factor = (numel(sp)/sum(smask(:)));
unmix = ismrm_calculate_sense_unmixing(acc_factor,smaps);

alias_img = ismrm_transform_kspace_to_image(data .* repmat(smask,[1 1 ncoils]),[1,2]);

%%
% Let's look at the aliased images:

ismrm_imshow(abs(alias_img),[],[2 4]);

%%
% It is evident already from this view that a couple of channels have a
% different noise level (channels 3 and 6). 

%%
% An now the actual SENSE reconstruction:

img = sum(unmix .* alias_img,3);
ismrm_imshow(abs(img))

%%
% In that reconstruction we did not take the noise correlation into
% account. Let's have a look at it. 

noise = reshape(noise_color,numel(noise_color)/ncoils, ncoils);
noise = permute(noise,[2 1]);
M = size(noise,2);
Rn = (1/(M-1))*(noise*noise');

figure;imagesc(abs(Rn)); axis equal; axis off; colormap(jet);

%%
% As we predicted from the raw aliased images, there is some increased
% noise in some channels. 
%
% Let's decoorelate:

dmtx = inv(chol(Rn,'lower'))*sqrt(2); %sqrt(2) for sd=1 in real and imag

%Both data and noise coil sensitivities
data_prew = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
smaps_prew = ismrm_apply_noise_decorrelation_mtx(smaps,dmtx);

%%
%Now let's repeat the reconstruction:
unmix_prew = ismrm_calculate_sense_unmixing(acc_factor,smaps_prew);

alias_img_prew = ismrm_transform_kspace_to_image(data_prew .* repmat(smask,[1 1 ncoils]),[1,2]);

img_prew = sum(unmix_prew .* alias_img_prew,3);
ismrm_imshow(cat(3,abs(img),abs(img_prew)));


%% SNR Scaled Reconstruction
%
% Now that we have performed noise-prewhitening, i.e. noise is scaled to
% sd=1, we would like to produce an image in SNR units. 
%
% For SNR scaled reconstruction, we need to make sure that all signal
% processing steps maintain SNR scaling. 
%
% Let's start by examining the transformation from k-space to image space:

%Generate white noise (like we have after prewhitening)
noise_white = complex(randn(size(data)),randn(size(data)));

%Mask out values that we don't acquire
noise_white = noise_white.*repmat(smask,[1 1 ncoils]);

%Transform from k-space to image:
noise_test = ismrm_transform_kspace_to_image(noise_white,[1,2]);

%%
% Get the standard deviation, we add a large value to simulate a high SNR
% situation.
% 
sd = std(abs(noise_test(:) + 1000))

%%
% As we can see the noise level is now about half of what it is supposed to
% be. This is because we have forgotten that only a quater of the samples
% are actually sampled and so we need to scale by the square root of the
% acceleration factor:

%Transform from k-space to image:
noise_test = sqrt(acc_factor)*ismrm_transform_kspace_to_image(noise_white,[1,2]);

%Now the standard deviation
sd = std(abs(noise_test(:) + 1000))

%%
% By scaling by the acceleration factor we have maintained unit noise
% scaling through our reconstruction. Now we can "trivially" obtain an SNR
% scaled reconstruction:

unmix_prew = ismrm_calculate_sense_unmixing(acc_factor,smaps_prew);

alias_img_prew = sqrt(acc_factor).*ismrm_transform_kspace_to_image(data_prew .* repmat(smask,[1 1 ncoils]),[1,2]);
img_snr = sum(unmix_prew .* alias_img_prew,3) ./ sqrt(sum(abs(unmix_prew).^2,3));

ismrm_imshow(abs(img_snr)); colorbar;

%% 
% We can also easily inspect the g-map directly
gmap = sqrt(sum(abs(unmix_prew).^2,3)).*sqrt(sum(abs(smaps_prew).^2,3));
ismrm_imshow(abs(gmap)); colormap(jet); colorbar;


%% Pseudo Replica Method
%
% What if we didn't have access to the unmixing coefficients directly?
% We can obtain SNR scaled reconstructions using the pseudo replica method.
% 

for r=1:100,
    noise_white = complex(randn(size(data)),randn(size(data)));
    noise_white = noise_white.*repmat(smask,[1 1 ncoils]);
    s = data_prew + noise_white;
    tmp = sqrt(acc_factor).*ismrm_transform_kspace_to_image(s .* repmat(smask,[1 1 ncoils]),[1,2]);
    img_noise_rep(:,:,r) = sum(tmp .* unmix_prew,3);
end

g_pseudo = std(abs(img_noise_rep + max(abs(img_noise_rep(:)))),[],3); 
g_pseudo(g_pseudo < eps) = 1;
snr_pseudo = mean(img_noise_rep,3)./g_pseudo;
g_pseudo = g_pseudo.*sqrt(sum(abs(smaps_prew).^2,3));
ismrm_imshow([abs(img_snr) abs(snr_pseudo)]); colorbar;
ismrm_imshow([abs(gmap) abs(g_pseudo)]); colormap(jet); colorbar;


%% Iterative Non-Cartesian SENSE