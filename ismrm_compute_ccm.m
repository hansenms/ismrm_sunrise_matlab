function ccm = ismrm_compute_ccm(csm, noise_matrix)
%
%  ccm = ismrm_compute_ccm(csm, noise_matrix)
%
%  Computes noise-optimal channel combination maps from  coil sensitivity
%  maps and a noise covariance matrix.
%
%  The ccm can be applied to channel-by-channel images as
%
%  im_composite = sum(ccm .* im_channel_by_channel, 3);
%
%  INPUT:
%    - csm [x,y, coil]          : coil sensitivity maps
%    - noise_matrix [coil,coil] : noise covariance matrix
%
%  OUTPUT:
%    - ccm [x,y, coil]          : channel combination maps
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%
nx = size(csm, 1);
ny = size(csm, 2);
nc = size(csm, 3);

if( nargin <2 || isempty(noise_matrix) )
    noise_matrix = eye(nc);
end

csm_matrix = reshape(csm, [nx*ny nc]);

relative_ccm = conj(csm_matrix) * pinv(noise_matrix);

scale_correction = abs(sum(relative_ccm .* csm_matrix, 2));

nonzero_ind = scale_correction>0;

ccm = zeros(size(csm_matrix));
ccm(nonzero_ind, :) = relative_ccm(nonzero_ind,:) ./ repmat(scale_correction(nonzero_ind), [1 nc]);

ccm = reshape(ccm, [nx ny nc]);