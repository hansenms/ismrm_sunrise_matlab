%%
% Demo of reconstruction in SNR units

%%
%Load Data
close all;
clear all;

load im1.mat
load smaps_phantom.mat

%Some settings
acc_factor = 4;
noise_level = 0.05*max(im1(:));
pseudo_replicas = 256;


%%
%Sample some data
noise_level = 0.05*max(im1(:));

[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 32);

noise = noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);
data_noise = data + noise;

noise = reshape(noise,size(noise,1)*size(noise,2),size(noise,3));
dmtx = ismrm_calculate_noise_decorrelation_mtx(noise(sp>0,:));

%Now decorrelation and after that sd is 1.0
data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
data_noise = ismrm_apply_noise_decorrelation_mtx(data_noise,dmtx);
smaps_prew = ismrm_apply_noise_decorrelation_mtx(smaps,dmtx);

[img_sense,g_sense,snr_sense,snr_pseudo_sense,g_pseudo_sense] = ismrm_cartesian_SENSE(data_noise .* repmat((sp == 1 | sp == 3),[1 1 size(smaps_prew,3)]),smaps_prew,acc_factor,pseudo_replicas);
[img_grappa,g_grappa,snr_grappa,snr_pseudo_grappa,g_pseudo_grappa] = ismrm_cartesian_GRAPPA(data_noise,sp,acc_factor,smaps_prew,pseudo_replicas);

ismrm_imshow([abs(img_sense) abs(img_grappa)]);
ismrm_imshow([abs(g_sense) abs(g_pseudo_sense) abs(g_grappa) abs(g_pseudo_grappa)]); colorbar;
ismrm_imshow([abs(snr_sense) abs(snr_pseudo_sense) abs(snr_grappa) abs(snr_pseudo_grappa)]); colorbar;