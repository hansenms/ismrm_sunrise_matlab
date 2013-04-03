%%
%Clean up
close all;
clear all

%%
%Load Data
load ../im1.mat
load ../smaps_phantom.mat
ncoils = size(smaps,3);
imsize = size(im1);


%%
%Show example data
close all;
figure(1);
imagescn(abs(im1));axis equal;axis off; colormap(gray);
figure(2)
imagescn(abs(smaps),[],[2 4])
figure(3)
imagescn(abs(smaps.*repmat(im1,[1 1 ncoils])),[],[2 4])

%%
%Undersampled image examples
close all;
data = smaps.*repmat(im1,[1 1 ncoils]);
ksp = ismrm_transform_image_to_kspace(data,[1,2]);
figure;imagescn(log(abs(ksp)),[],[2 4])
[data, sp] = ismrm_sample_data(permute(im1,[2 1]), permute(smaps,[2 1 3]), 2, 0);
figure;imagescn(abs(ismrm_transform_kspace_to_image(permute(data,[2 1 3]),[1,2])),[],[2 4]);
figure;imagescn(sqrt(sum(abs(ismrm_transform_kspace_to_image(permute(data,[2 1 3]),[1,2])).^2,3)));






%%
%Direct examplex
%

N = 64;

%Low res dataset
im1k = ismrm_transform_image_to_kspace(im1);
smapsk = ismrm_transform_image_to_kspace(smaps,[1,2]);

ra = (1:N)+(size(im1,1)-N)/2;

rho = ismrm_transform_kspace_to_image(im1k(ra, ra));
csm = ismrm_transform_kspace_to_image(smapsk(ra, ra,:),[1,2]);

tic;
fprintf('Forming encoding matrix....');

%Cartesian Fourier encoding
dftmtx2 = kron(dftmtx(N),dftmtx(N));

%Adding coil encoding 
for c=1:ncoils, 
    tmpc{c} = dftmtx2*diag(vec(csm(:,:,c))); 
    dftmtx2c = vertcat(tmpc{:}); 
end; clear tmpc;

%Acceleration factor 4
E = dftmtx2c(1:4:end,:);
%E = dftmtx2c;

fprintf('done\n');
toc;

tic;
fprintf('Sampling....');
s = E*rho(:);
fprintf('done\n');
toc;

%noise_level = 0.0*sqrt(sum(data.*conj(data)/numel(data)));
%s = s + noise_level*complex(randn(size(data)),randn(size(data)));

tic;
fprintf('Reconstruction....');
img = E \ s;
fprintf('done\n');
toc;
img = reshape(img,N,N);
showimage(img);


%%
%Simple GRAPPA and SENSE
close all;
acc_factor = 4;
noise_level = 0.05*max(im1(:));
[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 32);
data_noise = data + noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);

data       = data * acc_factor;
data_noise = data_noise .* acc_factor;

img_alias = ismrm_transform_kspace_to_image(data .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]),[1,2]);
img_alias_noise = ismrm_transform_kspace_to_image(data_noise .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]),[1,2]);
figure;showimage(img_alias);
figure;showimage(img_alias_noise);

[unmix_sense, gmap_sense]   = ismrm_calculate_sense_unmixing(acc_factor, smaps);
[unmix_grappa, gmap_grappa] = ismrm_calculate_grappa_unmixing(data_noise, [4 5], acc_factor, (sp > 1));

figure;
showimage(sum(img_alias .* unmix_sense,3),[2 3 1]); colorbar;title('SENSE (noiseless)');
showimage(sum(img_alias_noise .* unmix_sense,3),[2 3 2]); colorbar; title('SENSE (noise)'); 
showimage(gmap_sense,[2 3 3]); colorbar; title('SENSE g-factor');

showimage(sum(img_alias .* unmix_grappa,3),[2 3 4]); colorbar; title('GRAPPA (noiseless)');
showimage(sum(img_alias_noise .* unmix_grappa,3),[2 3 5]); colorbar; title('GRAPPA (noise)');
showimage(gmap_grappa,[2 3 6]); colorbar; title('GRAPPA g-factor');
colormap(gray);
set(gcf,'color','w');

figure;
showimage(sum(img_alias_noise .* unmix_sense,3),[2 2 1]); colorbar; title('SENSE'); 
showimage(gmap_sense,[2 2 2]); colorbar; title('SENSE g-factor');

showimage(sum(img_alias_noise .* unmix_grappa,3),[2 2 3]); colorbar; title('GRAPPA');
showimage(gmap_grappa,[2 2 4]); colorbar; title('GRAPPA g-factor');
colormap(gray);
set(gcf,'color','w');


%%
% SENSE example figures
close all;
acc_factor = 2;
[data, sp] = ismrm_sample_data(im1, smaps, acc_factor);
[data_full, sp_full] = ismrm_sample_data(im1, smaps, 1);
img_alias = ismrm_transform_kspace_to_image(data .* repmat(sp > 0,[1 1 size(smaps,3)]), [1,2]);


figure;
showimage(im1,[2 3 1]);axis off;
showimage(smaps(:,:,1),[2 3 2]);axis off;
showimage(ismrm_transform_kspace_to_image(data_full(:,:,1)),[2 3 3]);axis off;
showimage(log(abs(data_full(:,:,1))),[2 3 4]);axis off;
showimage(log(abs(data(:,:,1))),[2 3 5]);axis off;
showimage(ismrm_transform_kspace_to_image(data(:,:,1)),[2 3 6]);axis off;
set(gcf,'color','w');
colormap(gray);

figure;
showimage(im1,[2 3 1]);axis off;
showimage(smaps(:,:,1),[2 3 2]);axis off;
showimage(ismrm_transform_kspace_to_image(data(:,:,1)),[2 3 3]);axis off;

showimage(im1,[2 3 4]);axis off;
showimage(smaps(:,:,4),[2 3 5]);axis off;
showimage(ismrm_transform_kspace_to_image(data(:,:,4)),[2 3 6]);axis off;
set(gcf,'color','w');
colormap(gray);

%%
% Noise covariance
% Generate 8 channel examples
load ../noise_matrix_32ch.mat
Rn_normal = Rn;
load ../noise_matrix_32ch_broken.mat
Rn_broken = Rn;

figure;
showimage(Rn_normal / Rn_normal(1),[1 2 1]);caxis([0 2]);colorbar;
showimage(Rn_broken / Rn_normal(1),[1 2 2]);caxis([0 2]);colorbar;
set(gcf,'color','w');
colormap(jet);

L_normal = chol(Rn_normal,'lower');
L_broken = chol(Rn_broken,'lower');

%These are 32 channels, we should generate some 8 channel ones.
noise_samples = 10000;
noise_white_32 = complex(randn(size(Rn_normal,1),noise_samples),randn(size(Rn_normal,1),noise_samples));
noise_color_normal_32 = L_normal*noise_white_32;
noise_color_broken_32 = L_broken*noise_white_32;

noise_color_normal_8 = noise_color_normal_32(end-7:end,:);
noise_color_broken_8 = noise_color_broken_32(end-7:end,:);

Rn_normal_8 = (1/(noise_samples-1))*(noise_color_normal_8 * noise_color_normal_8');
Rn_normal_8 = Rn_normal_8 ./ Rn_normal_8(1);
Rn_broken_8 = (1/(noise_samples-1))*(noise_color_broken_8 * noise_color_broken_8');
Rn_broken_8 = Rn_broken_8 ./ Rn_broken_8(1);
L_normal_8 = chol(Rn_normal_8,'lower');
L_broken_8 = chol(Rn_broken_8,'lower');

figure;
showimage(Rn_broken_8,[1 2 2]);title('Normal coil')
caxis([0 3]);
showimage(Rn_normal_8,[1 2 1]);;title('Broken coil')
caxis([0 3]);
set(gcf,'color','w');
colormap(jet);


%%
%Simple SENSE with coloured noise
acc_factor = 4;
noise_level = 0.05*max(im1(:));
ncoils = size(smaps,3);
[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 32);
noise_white = noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);

noise_mtx{1} = zeros(size(smaps,3)); % No noise
noise_mtx{2} = eye(size(smaps,3));   % white noise
noise_mtx{3} = L_normal_8;           % Coil with some correlation
noise_mtx{4} = L_broken_8;           % Coil with some correlation

figure;
for n=1:length(noise_mtx),
    noise_color = reshape(permute(noise_mtx{n} * permute(reshape(noise_white, numel(data)/ncoils,ncoils),[2 1]),[2 1]),size(data));

    data_noise = data + noise_color;
    
    if (n > 1),
        data_prew  = reshape(permute(inv(noise_mtx{n}) * permute(reshape(data_noise, numel(data)/ncoils,ncoils),[2 1]),[2 1]),size(data));
        smaps_prew  = reshape(permute(inv(noise_mtx{n}) * permute(reshape(smaps, numel(data)/ncoils,ncoils),[2 1]),[2 1]),size(data));
    else
        data_prew   = data_noise;
        smaps_prew  = smaps;
    end        
    data_prew       = data_prew * acc_factor;
    data_noise      = data_noise .* acc_factor;

    img_alias_prew = ismrm_transform_kspace_to_image(data_prew .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]),[1,2]);
    img_alias_noise = ismrm_transform_kspace_to_image(data_noise .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]),[1,2]);

    [unmix_sense, gmap_sense]   = ismrm_calculate_sense_unmixing(acc_factor, smaps);
    [unmix_sense_prew, gmap_sense_prew]   = ismrm_calculate_sense_unmixing(acc_factor, smaps_prew);

    showimage(sum(img_alias_noise .* unmix_sense,3),[2, length(noise_mtx), n]);axis off;
    showimage(sum(img_alias_prew .* unmix_sense_prew,3),[2, length(noise_mtx), length(noise_mtx)+n]);axis off;
end

colormap(gray);
set(gcf,'color','w');

%%
%Recon in SNR Units
acc_factor = 4;
noise_level = 0.10*max(im1(:));
ncoils = size(smaps,3);
[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 32);
noise_white = noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);
noise_color = reshape(permute(L_broken_8 * permute(reshape(noise_white, numel(data)/ncoils,ncoils),[2 1]),[2 1]),size(data));

data_noise = data + noise_color;

eta = reshape(noise_color,numel(data)/ncoils,ncoils);
eta = permute(eta(find(sp > 0),:),[2 1]);

M = size(eta,2);
Psi = (1/(M-1))*(eta*eta');
L_inv = inv(chol(Psi,'lower'));

data_prew  = reshape(permute(L_inv * permute(reshape(data_noise, numel(data)/ncoils,ncoils),[2 1]),[2 1]),size(data));
smaps_prew  = reshape(permute(L_inv * permute(reshape(smaps, numel(data)/ncoils,ncoils),[2 1]),[2 1]),size(data));

data_prew       = data_prew * acc_factor;

img_alias_prew = ismrm_transform_kspace_to_image(data_prew .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]),[1,2]);
img_alias_noise = ismrm_transform_kspace_to_image(data_noise .* repmat(sp == 1 | sp == 3,[1 1 size(smaps,3)]),[1,2]);

[unmix, gmap]   = ismrm_calculate_sense_unmixing(acc_factor, smaps_prew);
%[unmix, gmap] = ismrm_calculate_grappa_unmixing(data_noise, [4 5], acc_factor, (sp > 1));

recon = sum(img_alias_prew .* unmix,3);
showimage(recon,[1 3 1]);colorbar;axis off;
showimage(gmap,[1 3 2]);colorbar;axis off;
dv = rssq(unmix,3);dv(dv < 1) = 1;
showimage(recon ./ dv,[1 3 3]);colorbar;axis off;
colormap(gray);
set(gcf,'color','w');

figure;
showimage(recon ./ dv);colorbar;axis off;
colormap(gray);
set(gcf,'color','w');

%%
%Cartesian SENSE with LSQR
close all;
acc_factor = 4;
noise_level = 0.05*max(im1(:));

if 0,
    load ../smaps_phantom.mat
    [data, sp] = ismrm_sample_data(im1, smaps, 1);
    noise = noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);
    data_noise = data + noise;
    coil_images = ismrm_transform_kspace_to_image(data_noise);
    smaps = ismrm_estimate_csm_walsh(coil_images);
    noise_level = 0.01*max(im1(:));
else
    load ../smaps_phantom.mat    
end

[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 32);

noise = noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);
data_noise = data + noise;

noise = reshape(noise,size(noise,1)*size(noise,2),size(noise,3));
dmtx = ismrm_calculate_noise_decorrelation_mtx(noise(sp>0,:));

data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
data_noise = ismrm_apply_noise_decorrelation_mtx(data_noise,dmtx);
smaps_prew = ismrm_apply_noise_decorrelation_mtx(smaps,dmtx);
csm_sq = sum(smaps_prew .* conj(smaps_prew),3); csm_sq(csm_sq < eps) = 1;

data       = data;
data_noise = data_noise;

samp_mat = (sp>0);%(sp == 1 | sp == 3);
s = data_noise(repmat(samp_mat,[1 1 size(smaps_prew,3)]) > 0);
E = @(x,tr) ismrm_encoding_cartesian_SENSE(x,smaps_prew,samp_mat,tr);
img_noise = lsqr(E, s, 1e-5,50);
img_noise = reshape(img_noise,size(smaps,1),size(smaps,2));
showimage(img_noise.*sqrt(csm_sq),[1 3 2]);colorbar;axis off;


%Pseudo-replica
if 0,
    image_formation_func = @(x) reshape(lsqr(E,x,1e-5,50),size(im1));
    s = data_noise(repmat(sp,[1 1 size(smaps_prew,3)]) > 0);
    reps = 25;
    [snr,g,noise_psf] = ismrm_pseudo_replica(s, image_formation_func,reps);
    img_noise_rep = reshape(img_noise_rep,[size(im1,1),size(im1,2), reps]);
    csm_sq = sum(smaps_prew .* conj(smaps_prew),3); csm_sq(csm_sq < eps) = 1;
    g = g .* sqrt(csm_sq);
    figure;
    showimage(snr,[1 2 1]); colorbar; axis off;
    showimage(g,[1 2 2]);colorbar; axis off;
end

[kx_cal,ky_cal] = ind2sub(size(sp),[find(sp > 1,1,'first') find(sp > 1,1,'last')]);

%Estimate SPIRiT kernel
kernel_mask = ones(7,7);
kernel_mask(4,4) = 0;

kernel = ismrm_estimate_convolution_kernel(data_noise(kx_cal(1):kx_cal(2),ky_cal(1):ky_cal(2),:), kernel_mask);
%kernel = ismrm_estimate_convolution_kernel(data_full, kernel_mask);

% nCoil = 8;
% CalibTyk = 0.02; 
% [AtA,] = corrMatrix(data_full,[7 7]);
% kernel = zeros([7 7 8 8]);
% for n=1:nCoil
%     disp(sprintf('Calibrating coil %d',n));
% 	kernel(:,:,:,n) = calibrate(AtA,[7 7],8,n,CalibTyk);
% end

padded_kernel = ismrm_transform_kernel_to_image_space(kernel,[size(im1,1), size(im1,2)]);


E = @(x,tr) ismrm_system_cartesian_SPIRiT(x,samp_mat,padded_kernel,tr);

img_spirit_noise = lsqr(E, [s;zeros(numel(smaps_prew),1)], 1e-5,50);
img_spirit_noise = reshape(img_spirit_noise,size(smaps_prew,1),size(smaps_prew,2),size(smaps_prew,3));
csm_sq = sum(smaps_prew .* conj(smaps_prew),3); csm_sq(csm_sq < eps) = 1;
img_spirit_noise = sum(conj(smaps_prew) .* img_spirit_noise,3) ./ csm_sq;
showimage(sqrt(csm_sq).*img_spirit_noise,[1 3 3]);colorbar;axis off;


img_alias_noise = sqrt(acc_factor)*ismrm_transform_kspace_to_image(data_noise .* repmat(sp == 1 | sp == 3,[1 1 size(smaps_prew,3)]),[1,2]);
[unmix_sense, gmap_sense]   = ismrm_calculate_sense_unmixing(acc_factor, smaps_prew);
showimage(sqrt(csm_sq).*sum(img_alias_noise .* unmix_sense,3),[1 3 1]);colorbar;axis off;

if 0,
    img_form_func = @(x) sum(ismrm_transform_kspace_to_image(sqrt(acc_factor)*x .* repmat(sp == 1 | sp == 3,[1 1 size(smaps_prew,3)]),[1,2]).*unmix_sense,3);
    [snr,g,noise_psf] = ismrm_pseudo_replica(data_noise, img_form_func,256);
    csm_sq = sum(smaps_prew .* conj(smaps_prew),3); csm_sq(csm_sq < eps) = 1;
    g = g .* sqrt(csm_sq);
    figure;
    showimage(snr,[1 3 1]); colorbar; axis off;
    showimage(g,[1 3 2]);colorbar; axis off;
    showimage(gmap_sense,[1 3 3]);colorbar; axis off;
end


colormap(gray);
set(gcf,'color','w');

%%
%Non-Cartesian (Radial) SENSE and SPIRiT
close all;
acc_factor = 8;

noise_level = 0.05*max(im1(:));

%Prepare NUFFT
N = [size(im1,1) size(im1,2)];
J = [5 5];
K = N*2;

[k,w] = ismrm_generate_radial_trajectory(size(im1,1), size(im1,1)*pi/2); %Fully sampled
area_weights = pi*(K(1)/2)^2; 
w = w .* (area_weights/sum(w(:)));

nufft_st = nufft_init(k*2*pi,N,J,K,N/2,'minmax:kb');

data_radial = nufft(repmat(im1,[1 1 size(smaps,3)]).*smaps,nufft_st)  ./ sqrt(prod(K));
data_radial = data_radial + noise_level*complex(randn(size(data_radial)),randn(size(data_radial)));

csm_sq = sum(smaps .* conj(smaps),3); csm_sq(csm_sq < eps) = 1;
recon_full = (sum(conj(smaps).*nufft_adj(data_radial .* repmat(w,[1 size(data_radial,2)]),nufft_st),3) ./ csm_sq) ./ sqrt(prod(K));

if 0,
    E = @(x,tr) ismrm_encoding_non_cartesian_SENSE(x,ones(size(im1)),nufft_st,w,tr);
    for c=1:size(smaps,3),
        cal_data(:,:,c) = reshape(lsqr(E, data_radial(:,c) .* sqrt(w), 1e-3,30),size(im1));
    end
    cal_data = ismrm_transform_image_to_kspace(cal_data,[1,2]);
end
%cal_data = ismrm_transform_image_to_kspace(nufft_adj(data_radial .* repmat(w,[1 size(data_radial,2)]),nufft_st) ./ sqrt(prod(K)),[1,2]);

%Now undersampled
projections = size(im1,1)/acc_factor;
[k,w] = ismrm_generate_radial_trajectory(size(im1,1), projections);

%Prepare NUFFT
N = [size(im1,1) size(im1,2)];
J = [5 5];
K = N*2;
nufft_st = nufft_init(k*2*pi,N,J,K,N/2,'minmax:kb');

%Scale weights such that the sum of the weights is the gridded area (oversampled)
area_weights = pi*(K(1)/2)^2; 
w = w .* (area_weights/sum(w(:)));

data_radial = nufft(repmat(im1,[1 1 size(smaps,3)]).*smaps,nufft_st)  ./ sqrt(prod(K));
data_radial = data_radial + noise_level*complex(randn(size(data_radial)),randn(size(data_radial)));

recon_undersampled = (1/sqrt(acc_factor))*(sum(conj(smaps).*nufft_adj(data_radial .* repmat(w,[1 size(data_radial,2)]),nufft_st),3) ./ csm_sq) ./ sqrt(prod(K));

E = @(x,tr) ismrm_encoding_non_cartesian_SENSE(x,smaps,nufft_st,w,tr);
img = lsqr(E, data_radial(:) .* repmat(sqrt(w),[size(smaps,3),1]), 1e-3,30);
img = reshape(img,size(smaps,1),size(smaps,2));


[kx_cal,ky_cal] = ind2sub(size(sp),[find(sp > 1,1,'first') find(sp > 1,1,'last')]);

%Estimate SPIRiT kernel
kernel_mask = ones(7,7);
kernel_mask(4,4) = 0;

%cal_data = ismrm_transform_image_to_kspace(nufft_adj(data_radial .* repmat(w,[1 size(data_radial,2)]),nufft_st) ./ sqrt(prod(K)),[1,2]);
%cal_data = data_noise;

cal_data = ismrm_transform_image_to_kspace(smaps,[1 2]);

cal_idx = 1:floor(projections/(pi/2));
cal_idx = cal_idx + bitshift(size(cal_data,1),-1)-bitshift(length(cal_idx),-1);

%kernel = ismrm_estimate_convolution_kernel(cal_data(cal_idx, cal_idx, :), kernel_mask);
kernel = ismrm_estimate_convolution_kernel(cal_data, kernel_mask);

padded_kernel = ismrm_transform_kernel_to_image_space(kernel,[size(im1,1), size(im1,2)]);

E = @(x,tr) ismrm_system_non_cartesian_SPIRiT(x,nufft_st,w(:),padded_kernel,tr);

img_spirit = lsqr(E, [data_radial(:) .* repmat(sqrt(w),[size(smaps,3),1]);zeros(numel(smaps),1)], 1e-5,50);
img_spirit = reshape(img_spirit,size(smaps,1),size(smaps,2),size(smaps,3));
csm_sq = sum(smaps .* conj(smaps),3); csm_sq(csm_sq < eps) = 1;
img_spirit = sum(conj(smaps) .* img_spirit,3) ./ csm_sq;

showimage(recon_full,[1 4 1]);colorbar;axis off;
showimage(recon_undersampled,[1 4 2]);colorbar;axis off;
showimage(img,[1 4 3]);colorbar;axis off;
showimage(img_spirit,[1 4 4]);colorbar;axis off;

colormap(gray);
set(gcf,'color','w');





