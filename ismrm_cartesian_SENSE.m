function [img,gmap,snr,snr_pseudo,gmap_pseudo,noise_psf_pseudo] = ismrm_cartesian_SENSE(inp,csm,acc_factor,replicas)
%
%   [img,gmap,snr_pseudo,gmap_pseudo,noise_psf_pseudo] = ismrm_cartesian_SENSE(inp,csm,acc_factor,replicas)
%
%   Cartesian SENSE reconstruction.
%
%   It is recommended to input data with pre-whitened noise scale to sd=1.
%
%
%   INPUT:
%     - inp         [kx,ky,coils]        : Undersampled input k-space data
%     - acc+factor  scalar               : Parallel Imaging Acceleration factor
%     - csm         [x,y,coil]           : Optional coil sensitivity map (for coil combination)
%     - replicas    scalar (dafault 100) : Number of replicas to run if SNR
%                                          is requested
%
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

if nargin<4,
    replicas = 100;
end

[unmix_sense, gmap]   = ismrm_calculate_sense_unmixing(acc_factor, csm);


img_alias = sqrt(acc_factor)*ismrm_transform_kspace_to_image(inp,[1,2]);
img = sum(img_alias .* unmix_sense,3);
snr = img ./ sqrt(sum(abs(unmix_sense).^2,3));


if (nargout > 3),
    image_formation_func = @(x) sum(ismrm_transform_kspace_to_image(sqrt(acc_factor)* inp,[1,2]).*unmix_sense,3);
    figure;showimage(image_formation_func(inp));colorbar;
    [snr_pseudo,gmap_pseudo,noise_psf_pseudo] = ismrm_pseudo_replica(inp, image_formation_func,replicas);
    
    gmap_pseudo = gmap_pseudo .* sqrt(sum(abs(csm).^2,3));
end

