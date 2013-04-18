function [csm] = ismrm_estimate_csm_mckenzie(img)
%
%   [csm] = ismrm_estimate_csm_mckenzie(img)
%
%   Estimates relative coil sensitivity maps from a set of
%   channel-by-channel images.
%
%   McKenzie et al. (Magn Reson Med
%   2002;47:529-538.)
%
%   INPUT:
%     - img     [x,y,coil]   : Coil images
%
%   OUTPUT:
%
%     - csm     [x,y,coil]    : Relative coil sensitivity maps
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%

nx = size(img,1);
ny = size(img,2);
nc = size(img,3);

img_matrix = reshape(img, [nx*ny nc]);
scale_correction = sqrt(sum(abs(img_matrix).^2,2));
nonzero_ind = scale_correction > 0;
csm = zeros(size(img_matrix));
csm(nonzero_ind,:) = img_matrix(nonzero_ind,:) ./ repmat( scale_correction(nonzero_ind), [1 nc]);
csm = reshape(csm, [nx ny nc]);

