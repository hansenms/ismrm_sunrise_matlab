function aem = ismrm_calculate_aem(pixel_mask, true_csm, unmix_im, acc_factor)
%
%  ismrm_calculate_aem(pixel_mask, true_csm, unmix_im, acc_factor)
%
%  Computes the square root of an "aliasing energy map" 
%
%  INPUT:
%    - pixel_mask [x,y]      : 1 = pixel that might have signal
%                              0 = pixel that won't have signal
%    - true_csm [x,y,coil]   : ground truth coil sensitivity maps
%    - unmix_im [x,y,coil]   : unmix images under evaluation
%    - acc_factor [int]      : acceleration factor corresponding to
%                              unmix_im
%
%  OUTPUT:
%    - aem [x,y]             : square root of aliasing energy map
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%

im_shape = size(pixel_mask);

ncoils = size(true_csm, 3);
masked_csm = true_csm .* repmat(pixel_mask, [1 1 ncoils]);

aem = zeros(im_shape);

for a_ind = 1:acc_factor,
    sampled_fov = im_shape(2)/acc_factor;
    signal_ind = (1:sampled_fov) + (a_ind-1)*sampled_fov;
    im_alias = repmat(masked_csm(:, signal_ind,:), [1 acc_factor, 1]);

    im_hat = abs(sum(im_alias .* unmix_im, 3)) ./ acc_factor;

%    im_start = zeros(im_shape);
%    im_start(:,signal_ind) = pixel_mask(:,signal_ind);
%    ismrm_imshow(im_start, [0 1]);
%    ismrm_imshow(im_hat, [0 .1]);

    im_hat(:,signal_ind) = 0;
    aem = aem + im_hat.^2;
end

aem = sqrt(aem);