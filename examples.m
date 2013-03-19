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
%Simple SENSE
close all;
[data, sp] = ismrm_sample_data(im1, smaps, 4, 0);
img_alias = ismrm_transform_kspace_to_image(data,[1,2]);
figure;showimage(img_alias);
unmix = zeros(size(im1,1),size(im1,2),ncoils); 
for x=1:size(im1,1), 
    unmix(x,:,:) = ismrm_calculate_sense_unmixing_1d(4, squeeze(smaps(x,:,:))); 
end
figure;showimage(sum(img_alias .* unmix,3));

