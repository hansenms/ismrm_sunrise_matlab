%%
% Simple demo of iterative parallel imaging using the iteratieve SENSE
% and SPIRiT methods

%%
%Load Data
close all;
clear all;

load im1.mat
load smaps_phantom.mat

%Some settings
acc_factor = 4;
noise_level = 0.05*max(im1(:));
pseudo_replicas = 0; %zero means no pseudo replicas will be done. 


%%
%Sample some data
noise_level = 0.05*max(im1(:));

[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 32);

noise = noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);
data_noise = data + noise;

noise = reshape(noise,size(noise,1)*size(noise,2),size(noise,3));
dmtx = ismrm_calculate_noise_decorrelation_mtx(noise(sp>0,:));

data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
data_noise = ismrm_apply_noise_decorrelation_mtx(data_noise,dmtx);
smaps_prew = ismrm_apply_noise_decorrelation_mtx(smaps,dmtx);

%%
%Now some reconstructions

%Create a blurred regularization image
if (0),
    reg_size = 32;
    reg_img = ismrm_transform_kspace_to_image(padarray(hamming(reg_size)*hamming(reg_size)',[(size(im1,1)-reg_size)/2 (size(im1,2)-reg_size)/2],0,'both').*ismrm_transform_image_to_kspace(im1));
    reg_img = abs(reg_img)+1;
else
    reg_img = [];
end

samp_mat = (sp == 1 | sp == 3);
s = data_noise(repmat(samp_mat,[1 1 size(smaps_prew,3)]));

if (pseudo_replicas > 0),
    [img_it_sense,snr_it_sense,g_it_sense] = ismrm_cartesian_iterative_SENSE(s,samp_mat,smaps_prew,reg_img,pseudo_replicas);
else
    [img_it_sense] = ismrm_cartesian_iterative_SENSE(s,samp_mat,smaps_prew,reg_img);
end

[kx_cal,ky_cal] = ind2sub(size(sp),[find(sp > 1,1,'first') find(sp > 1,1,'last')]);
cal_data = data_noise(kx_cal(1):kx_cal(2),ky_cal(1):ky_cal(2),:);

if (pseudo_replicas > 0),
    [img_spirit,snr_spirit,g_spirit] = ismrm_cartesian_SPIRiT(s,samp_mat,cal_data,smaps_prew,pseudo_replicas);
else
    [img_spirit] = ismrm_cartesian_SPIRiT(s,samp_mat,cal_data,smaps_prew);
end

if (pseudo_replicas > 0),
    [img_sense,g_sense,snr_sense,snr_pseudo_sense,g_pseudo_sense] = ismrm_cartesian_SENSE(data_noise .* repmat((sp == 1 | sp == 3),[1 1 size(smaps_prew,3)]),smaps_prew,acc_factor,pseudo_replicas);
    [img_grappa,g_grappa,snr_grappa,snr_pseudo_grappa,g_pseudo_grappa] = ismrm_cartesian_GRAPPA(data_noise,sp,acc_factor,smaps_prew,pseudo_replicas);
else
    [img_sense,g_sense,snr_sense] = ismrm_cartesian_SENSE(data_noise .* repmat((sp == 1 | sp == 3),[1 1 size(smaps_prew,3)]),smaps_prew,acc_factor,pseudo_replicas);
    [img_grappa,g_grappa,snr_grappa] = ismrm_cartesian_GRAPPA(data_noise,sp,acc_factor,smaps_prew);
end

ismrm_imshow([abs(img_sense) abs(img_grappa) abs(img_it_sense) abs(img_spirit)]);
if (pseudo_replicas > 0),
    ismrm_imshow([abs(g_sense) abs(g_grappa) abs(g_it_sense) abs(g_spirit)]);
else
    ismrm_imshow([abs(g_sense) abs(g_grappa)]);
end