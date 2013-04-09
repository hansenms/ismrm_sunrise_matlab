function [img,gmap,snr,snr_pseudo,gmap_pseudo,noise_psf_pseudo] = ismrm_cartesian_GRAPPA(inp,samp_mat, acc_factor,csm,replicas)
%
%   [img,gmap,snr_pseudo,gmap_pseudo,noise_psf_pseudo] = ismrm_cartesian_GRAPPA(inp,csm,acc_factor,replicas)
%
%   Cartesian SENSE reconstruction.
%
%   It is recommended to input data with pre-whitened noise scale to sd=1.
%
%
%   INPUT:
%     - inp         [kx,ky,coils]        : Undersampled input k-space data
%     - acc_factor  scalar               : Parallel Imaging Acceleration factor
%     - csm         [x,y,coil]           : Optional coil sensitivity map (for coil combination)
%     - replicas    scalar (dafault 100) : Number of replicas to run if SNR
%                                          is requested
%     - samp_mat     [kx,ky,c]           : Sampling pattern (0 = not sampled,
%                                           1 = imaging data,
%                                           2 = reference data,
%                                           3 = reference and imaging data)%
%   OUTPUT:
%     - img              [x,y]                : Output image
%     - gmap             [x,y]                : g-map (calulated from unmixing coefficients)
%     - snr              [x,y]                : An image in SNR units. From direct recon.
%     - snr_pseudo       [x,y]                : An image in SNR units. Using psudo replica method.
%     - gmap_pseudo      [x,y]                : A g-map (Pseudo replica)
%     - noise_psf_pseudo [x,y]                : Point spread function of the noise
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%

if nargin<5,
    replicas = 100;
end

if nargin<4,
   csm = []; 
end

[unmix, gmap] = ismrm_calculate_grappa_unmixing(inp, [5 4], acc_factor, (samp_mat > 1),csm);


img_alias = sqrt(acc_factor)*ismrm_transform_kspace_to_image(inp .* repmat((samp_mat == 1 | samp_mat == 3),[1 1 size(inp,3)]),[1,2]);
img = sum(img_alias .* unmix,3);
snr = img ./ sqrt(sum(abs(unmix).^2,3));


if (nargout > 3),
    image_formation_func = @(x) sum(ismrm_transform_kspace_to_image(sqrt(acc_factor)*x .* repmat((samp_mat == 1 | samp_mat == 3),[1 1 size(csm,3)]),[1,2]).*unmix,3);
    [snr_pseudo,gmap_pseudo,noise_psf_pseudo] = ismrm_pseudo_replica(inp, image_formation_func,replicas);
    gmap_pseudo = gmap_pseudo .* sqrt(sum(abs(csm).^2,3));
end

