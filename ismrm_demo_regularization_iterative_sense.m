%%
% Simple demo of non-Cartesian parallel imaging with regularization

%%
%Load Data
close all;
clear all;

load im1.mat
load smaps_phantom.mat

%Some settings
acc_factor = 4;
noise_level = 0.30*max(im1(:));
trajectory = 'spiral';
pseudo_replicas = 0; %zero means no pseudo replicas will be done. 
lambda = [1.2]; %Regularization
%lambda = [0.1 0.5 0.8 1.0 1.2 1.5 2.0 5.0]; %Regularization

%%
% Simlate data
%Non-Cartesian (Radial) SENSE and SPIRiT

if (strcmp(trajectory,'radial')),
    projections = floor(size(im1,1)/acc_factor);
    [k,w] = ismrm_generate_radial_trajectory(size(im1,1), projections);
    area_weights = pi*(0.5)^2; %Fraction of k-space covered
    w = w .* (area_weights/sum(w(:)));    
elseif (strcmp(trajectory,'random'))
    k = rand(numel(im1)/acc_factor,2)-0.5;
    w = ones(size(k,1),1);
    area_weights = 1;
    w = w .* (area_weights/sum(w(:)));        
elseif (strcmp(trajectory,'spiral')),
    [k,w] = ismrm_calculate_spiral_trajectory(acc_factor);
else
    error('Invalid trajectory')
end


%Prepare NUFFT
N = [size(im1,1) size(im1,2)];
J = [5 5];
K = N*2;
nufft_st = nufft_init(k*2*pi,N,J,K,N/2,'minmax:kb');

%Sample
data_noiseless = nufft(repmat(im1,[1 1 size(smaps,3)]).*smaps,nufft_st)  ./ sqrt(prod(N));
noise = noise_level*complex(randn(size(data_noiseless)),randn(size(data_noiseless)));
data = data_noiseless + noise;

%Prewhiten data, after this noise sd = 1.0
dmtx = ismrm_calculate_noise_decorrelation_mtx(noise);
data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
smaps_prew = ismrm_apply_noise_decorrelation_mtx(smaps,dmtx);


%Create a blurred regularization image
reg_size = 32;
reg_img = ismrm_transform_kspace_to_image(padarray(hamming(reg_size)*hamming(reg_size)',[(size(im1,1)-reg_size)/2 (size(im1,2)-reg_size)/2],0,'both').*ismrm_transform_image_to_kspace(im1));
reg_img = abs(reg_img)+1;

%%
% Reconstruct undersampled data

if (pseudo_replicas > 1),
    [img_reg,snr_reg,g_reg,noise_psf_reg] = ismrm_non_cartesian_sense(data,k,w,smaps_prew,reg_img,lambda,pseudo_replicas);
    [img_sense,snr_sense,g_sense,noise_psf_sense] = ismrm_non_cartesian_sense(data,k,w,smaps_prew,[],lambda,pseudo_replicas);
else
    for la=1:length(lambda),
        [img_reg(:,:,la)] = ismrm_non_cartesian_sense(data(:),k,w,smaps_prew,reg_img,lambda(la));
    end
    [img_sense] = ismrm_non_cartesian_sense(data(:),k,w,smaps_prew,[],lambda);
end

ismrm_imshow(cat(3,abs(img_sense), abs(img_reg(:,:,1)))); colormap(gray);
if (pseudo_replicas > 1),
    ismrm_imshow(cat(3,abs(snr_sense), abs(snr_reg))); colormap(gray);
    ismrm_imshow(cat(3,abs(g_sense), abs(g_reg))); colormap(jet);    
end

if (length(lambda) > 1),
    ismrm_imshow(abs(img_reg),[],[2 4]); colormap(gray);
end


