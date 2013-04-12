%%
% Demo of effects of noise correlation

%%
%Load Data
close all;
clear all;

load im1.mat
load smaps_phantom.mat
load noise_covariances.mat

L_normal_8 = chol(Rn_normal_8,'lower');
L_broken_8 = chol(Rn_broken_8,'lower');

%Some settings
acc_factor = 4;
noise_level = 0.05*max(im1(:));
pseudo_replicas = 256;


%%
%Sample some data%%
%Simple SENSE with coloured noise
ncoils = size(smaps,3);
[data, sp] = ismrm_sample_data(im1, smaps, acc_factor, 32);
noise_white = noise_level*complex(randn(size(data)),randn(size(data))) .* repmat(sp > 0,[1 1 size(smaps,3)]);

noise_mtx{1} = eye(size(smaps,3));   % white noise
noise_mtx{2} = L_normal_8;           % Coil with some correlation
noise_mtx{3} = L_broken_8;           % Coil with some correlation

for n=1:length(noise_mtx),
    noise_color = reshape(permute(noise_mtx{n} * permute(reshape(noise_white, numel(data)/ncoils,ncoils),[2 1]),[2 1]),size(data));

    data_noise = data + noise_color;
    
    noise = reshape(noise_color,size(noise_color,1)*size(noise_color,2),size(noise_color,3));
    dmtx = ismrm_calculate_noise_decorrelation_mtx(noise(sp>0,:));

    %Now decorrelation and after that sd is 1.0
    data_prew = ismrm_apply_noise_decorrelation_mtx(data_noise,dmtx);
    smaps_prew = ismrm_apply_noise_decorrelation_mtx(smaps,dmtx);

    [img_sense_prew{n},g_sense_prew{n}] = ismrm_cartesian_SENSE(data_prew .* repmat((sp == 1 | sp == 3),[1 1 size(smaps_prew,3)]),smaps_prew,acc_factor);
    [img_grappa_prew{n},g_grappa_prew{n}] = ismrm_cartesian_GRAPPA(data_prew,sp,acc_factor,smaps_prew);

    [img_sense{n},g_sense{n}] = ismrm_cartesian_SENSE(data_noise .* repmat((sp == 1 | sp == 3),[1 1 size(smaps,3)]),smaps,acc_factor);
    [img_grappa{n},g_grappa{n}] = ismrm_cartesian_GRAPPA(data_noise,sp,acc_factor,smaps);
end


ismrm_imshow(cat(3,abs(img_sense{1}), abs(img_sense{2}), abs(img_sense{3}),abs(img_sense_prew{1}), abs(img_sense_prew{2}), abs(img_sense_prew{3})),[],[2 3]);
ismrm_imshow(cat(3,abs(img_grappa{1}), abs(img_grappa{2}), abs(img_grappa{3}),abs(img_grappa_prew{1}), abs(img_grappa_prew{2}), abs(img_grappa_prew{3})),[],[2 3]);

ismrm_imshow(cat(3,abs(g_sense{1}), abs(g_sense{2}), abs(g_sense{3}),abs(g_sense_prew{1}), abs(g_sense_prew{2}), abs(g_sense_prew{3})),[1 5],[2 3]);colormap(jet);
ismrm_imshow(cat(3,abs(g_grappa{1}), abs(g_grappa{2}), abs(g_grappa{3}),abs(g_grappa_prew{1}), abs(g_grappa_prew{2}), abs(g_grappa_prew{3})),[1 5],[2 3]);colormap(jet);
