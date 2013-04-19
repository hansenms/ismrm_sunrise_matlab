function [im_out, correction_image] = ismrm_normalize_shading_to_sos(im_in)
%
%  ismrm_normalize_shading_to_sos(im_in)
%
%  Applies correction to csm or ccm images so that the shading profile is
%  the same as a square root sum-of-squares channel combination.  This
%  allows normalization of the shading profile between different
%  reconstruction methods.
%
%  INPUT:
%    - im_in [x, y, coils] : input ccm or csm images
%
%  OUTPUT:
%    - im_out [x, y, coils] : shading normalized output ccm or csm
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%
nx = size(im_in, 1);
ny = size(im_in, 2);
nc = size(im_in, 3);

im_in_matrix = reshape(im_in, [nx*ny nc]);
shading_correction = sqrt(sum(abs(im_in_matrix).^2, 2));
nonzero_ind = shading_correction > 0;

correction_image = zeros(size(shading_correction));
correction_image(nonzero_ind) = 1 ./ shading_correction(nonzero_ind);
correction_image = reshape(correction_image, [nx ny]);

im_out = im_in .* repmat(correction_image, [1 1 nc]);
%im_out = zeros(size(im_in_matrix));
%im_out(nonzero_ind,:) = im_in_matrix(nonzero_ind,:) ./ repmat(shading_correction(nonzero_ind), [1 nc]);

%im_out = reshape(im_out, [nx ny nc]);