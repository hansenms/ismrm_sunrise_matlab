%%
%Load and sample data
clear all; close all;
load ../im1.mat
load ../smaps_phantom.mat

%Set basic simumation parameters
accel = 4;
ref_lines = 32;

noise_level = 0.05;

[data, pat] = ismrm_sample_data(im1,smaps,accel,ref_lines);
[data_full, pat_full] = ismrm_sample_data(im1,smaps,1,ref_lines);

[kx_cal,ky_cal] = ind2sub(size(pat),[find(pat > 1,1,'first') find(pat > 1,1,'last')]);

noise_magnitude = noise_level*max(abs(im1(:)));

data = data + noise_magnitude*complex(randn(size(data)),randn(size(data))).*repmat(pat>0,[1 1 size(data,3)]);


%Let's create a simple radial sampling pattern:
projections = 32;
kx = linspace(-pi,pi,size(data,1));
om = zeros(projections*numel(kx),2);
weights = zeros(size(om,1),1);
for p=1:projections,
    a = p*pi/(projections-1);
    om((p-1)*numel(kx)+(1:numel(kx)),1) = kx'*cos(a);
    om((p-1)*numel(kx)+(1:numel(kx)),2) = kx'*sin(a);
    weights((p-1)*numel(kx)+(1:numel(kx))) = abs(kx);
end
weights = weights .* (numel(weights)/sum(weights(:)));

N = [size(im1,1) size(im1,2)];
J = [5 5];
K = N*2;
nufft_st = nufft_init(om,N,J,K,N/2,'minmax:kb');
data_radial = nufft(repmat(im1,[1 1 size(smaps,3)]).*smaps,nufft_st)  ./ sqrt(prod(K));
data_radial = data_radial + noise_magnitude*complex(randn(size(data_radial)),randn(size(data_radial)));

csm_sq = sum(smaps .* conj(smaps),3) + eps;


%%
%Cartesian GRAPPA reconstruction

[unmix, csm] = ismrm_calculate_grappa_unmixing(data, [5 4], 4, (pat > 1));
recon_cartesian_grappa = sum(ismrm_transform_kspace_to_image(data .* repmat((pat == 1 | pat == 3),[1 1 size(data,3)]),[1,2]) .* unmix,3);
showimage(recon_cartesian_grappa);title('Cartesian GRAPPA');

%%
%Simple conjugate gradient sense
m = data(repmat((pat == 1 | pat == 3),[1 1 size(data,3)]));
E = @(x,trsn) ismrm_encoding_cartesian_SENSE(x,smaps,(pat == 1 | pat == 3),trsn);
[im,FLAG,RELRES,ITER,RESVEC] = lsqr(E,m,1e-3,50);
im = reshape(im,size(smaps,1),size(smaps,2));

im2 = ismrm_cg(m(:), E,'print_residual',1,'limit',1e-6,'iterations',50);
im2 = reshape(im2,size(smaps,1),size(smaps,2));

%%
%Estimate SPIRiT kernel
kernel_mask = ones(7,7);
kernel_mask(4,4) = 0;

kernel = ismrm_estimate_convolution_kernel(data(kx_cal(1):kx_cal(2),ky_cal(1):ky_cal(2),:), kernel_mask);
%kernel = ismrm_estimate_convolution_kernel(data_full, kernel_mask);

% nCoil = 8;
% CalibTyk = 0.02; 
% [AtA,] = corrMatrix(data_full,[7 7]);
% kernel = zeros([7 7 8 8]);
% for n=1:nCoil
%     disp(sprintf('Calibrating coil %d',n));
% 	kernel(:,:,:,n) = calibrate(AtA,[7 7],8,n,CalibTyk);
% end


kernel = flipdim(flipdim(kernel,1),2);

% data = zeros([N 5]);
padded_kernel = zeros(size(data,1),size(data,2),size(data,3),size(data,3));
padded_kernel([1:size(kernel,1)]+bitshift(size(data,1)-size(kernel,1)-1,-1)+1, ...
    [1:size(kernel,2)]+bitshift(size(data,2)-size(kernel,2)-1,-1)+1, :, :) = kernel;

padded_kernel = fftshift(ifft(ifftshift(padded_kernel,1),[],1),1);
padded_kernel = fftshift(ifft(ifftshift(padded_kernel,2),[],2),2);
padded_kernel = padded_kernel*(size(padded_kernel,1)*size(padded_kernel,2));


L = @(x,trsn) ismrm_image_kernel_SPIRiT(x,padded_kernel,trsn);
E = @(x,trsn) ismrm_encoding_cartesian_SPIRiT(x,(pat == 1 | pat == 3),trsn);
im = ismrm_cg(m(:), E,'print_residual',1,'limit',1e-5,'fL',L,'iterations',50);
im = reshape(im,size(smaps,1),size(smaps,2),size(smaps,3));
im = sum((im.*conj(smaps)),3) ./ csm_sq;
showimage(im)

A = @(x,trsn) ismrm_system_cartesian_SPIRiT(x,(pat == 1 | pat == 3),padded_kernel,trsn);
imgelements = numel(smaps);
[im,FLAG,RELRES,ITER,RESVEC] = lsqr(A,[m;zeros(imgelements,1)],1e-3,20);
im = reshape(im,size(smaps,1),size(smaps,2),size(smaps,3));
im = sum((im.*conj(smaps)),3) ./ csm_sq;
showimage(im)

%%
%Non-Cartesian SPIRiT
E = @(x,trsn) ismrm_sampling_non_cartesian_SPIRiT(x,nufft_st,weights,trsn);
%im = ismrm_cg(data_radial(:), E,'print_residual',1,'limit',1e-5,'iterations',50,'fL',L);
im = ismrm_cg(vec(data_radial.*repmat(sqrt(weights),[1 size(smaps,3)])), E,'print_residual',1,'limit',1e-5,'iterations',50,'fL',L);
im = reshape(im,size(smaps,1),size(smaps,2),size(smaps,3));
im = sum((im.*conj(smaps)),3) ./ csm_sq;
showimage(im,[1 2 1]);title('Radial SPIRiT');

%E = @(x,trsn) ismrm_sampling_non_cartesian_SPIRiT(x,nufft_st,weights,trsn);
%im = ismrm_cg(vec(data_radial.*repmat(sqrt(weights),[1 size(smaps,3)])), E,'print_residual',1,'limit',1e-5,'iterations',50,'fL',L);
A = @(x,trsn) ismrm_system_non_cartesian_SPIRiT(x,nufft_st,weights,padded_kernel,trsn);
imgelements = numel(smaps);
[im,FLAG,RELRES,ITER,RESVEC] = lsqr(A,[vec(data_radial.*repmat(sqrt(weights),[1 size(smaps,3)]));zeros(imgelements,1)],1e-3,20);
im = reshape(im,size(smaps,1),size(smaps,2),size(smaps,3));
im = sum((im.*conj(smaps)),3) ./ csm_sq;
showimage(im,[1 2 1]);title('Radial SPIRiT');


%Now with non-Cartesian SENSE
E = @(x,trsn) ismrm_encoding_non_cartesian_SENSE(x,smaps,nufft_st,weights,trsn);
%im  = ismrm_cg(data_radial(:), E,'print_residual',1,'limit',1e-3,'iterations',20,'lambda',1);
im  = ismrm_cg(vec(data_radial.*repmat(sqrt(weights),[1 size(smaps,3)])), E,'print_residual',1,'limit',1e-3,'iterations',20,'lambda',1);

im = reshape(im,size(smaps,1),size(smaps,2));
showimage(im,[1 2 2]);title('Radial CG SENSE');


