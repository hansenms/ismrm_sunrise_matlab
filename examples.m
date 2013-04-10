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
% Noise covariance
% Generate 8 channel examples
load ../noise_matrix_32ch.mat
Rn_normal = Rn;
load ../noise_matrix_32ch_broken.mat
Rn_broken = Rn;

figure;
showimage(Rn_normal / Rn_normal(1),[1 2 1]);caxis([0 2]);colorbar;
showimage(Rn_broken / Rn_normal(1),[1 2 2]);caxis([0 2]);colorbar;

figure;
showimage(Rn_broken / Rn_normal(1),[1 2 1]);caxis([0 2]);colorbar;
showimage(eye(32),[1 2 2]);caxis([0 2]);colorbar;


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









