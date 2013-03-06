clear all; close all;
load ../im1.mat
load ../smaps_phantom.mat


accel = 4;
ref_lines = 32;
noise_level = 0.05;


[data, pat] = ismrm_sample_data(im1,smaps,accel,ref_lines);
[data_full, pat_full] = ismrm_sample_data(im1,smaps,1,ref_lines);

[unmix, csm] = ismrm_calculate_grappa_unmixing(data, [5 4], 4, (pat_full > 1));
recon = sum(ismrm_transform_kspace_to_image(data .* repmat((pat == 1 | pat == 3),[1 1 size(data,3)]),[1,2]) .* unmix,3);


[kx_cal,ky_cal] = ind2sub(size(pat),[find(pat > 1,1,'first') find(pat > 1,1,'last')]);


%Simple conjugate gradient sense
m = data(repmat((pat == 1 | pat == 3),[1 1 size(data,3)]));
E = @(x,trsn) ismrm_encoding_cartesian_SENSE(x,smaps,(pat == 1 | pat == 3),trsn);
im = lsqr(E,m,1e-3,20);
im = reshape(im,size(csm,1),size(csm,2));

%Estimate SPIRiT kernel
kernel_mask = ones(7,7);
kernel_mask(4,4) = 0;

kernel = ismrm_estimate_convolution_kernel(data(kx_cal(1):kx_cal(2),ky_cal(1):ky_cal(2),:), kernel_mask);
kernel = flipdim(flipdim(kernel,1),2);
padded_kernel = zeros(size(data,1),size(data,2),size(data,3),size(data,3));
padded_kernel([1:size(kernel,1)]+bitshift(size(data,1)-size(kernel,1)-1,-1)+1, ...
    [1:size(kernel,2)]+bitshift(size(data,2)-size(kernel,2)-1,-1)+1, :, :) = kernel;

padded_kernel = fftshift(ifft(ifftshift(padded_kernel,1),[],1),1);
padded_kernel = fftshift(ifft(ifftshift(padded_kernel,2),[],2),2);
padded_kernel = padded_kernel*(size(padded_kernel,1)*size(padded_kernel,2));


L = @(x,trsn) ismrm_image_kernel_SPIRiT(x,padded_kernel,trsn);
E = @(x,trsn) ismrm_encoding_cartesian_SPIRiT(x,(pat == 1 | pat == 3),trsn);
[rho] = ismrm_cg(m(:), E,'print_residual',1,'limit',1e-5,'fL',L,'iterations',50);




