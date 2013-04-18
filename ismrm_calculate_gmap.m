function gmap = ismrm_calculate_gmap(unmixing, ccm, noise_matrix, acc_factor)
%
%  ismrm_calculate_gmap(unmixing, ccm, noise_matrix, acc_factor)
%
%  Computes g-factor map (relative noise enhancement between unaccelerated
%  and accelerated case normalized by scan time.
%
%  N.B. This function assumes that unmixing operates on aliased images
%  formed from data that is scaled the same as the fully sampled case, but
%  with (R-1)/R of the data set to zero.  As such, std dev. noise of fully
%  sampled channel-by-channel images = sqrt(R) * std dev. noise of
%  channel-by-channel aliased images.  This is accounted for by dividing by
%  R in this function instead of sqrt(R).  When the unmixing and ccm are
%  applied to the accelerated and fully sampled data, the final images
%  should be scaled identically.
%
%  INPUT:
%    - unmixing [x,y, coil]     : unmixing images for accelerated case.
%    - ccm [x,y,coil]           : channel combination maps
%    - noise_matrix [coil,coil] : noise covariance matrix
%    - acc_factor [int]      : acceleration factor corresponding to
%                              unmixing
%
%  OUTPUT:
%    - gmap [x,y]             : g-factor map
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%
if nargin < 3,
    noise_matrix = [];
end
if isempty(noise_matrix),
    noise_matrix = eye(size(unmixing,3));
end

nx = size(unmixing,1);
ny = size(unmixing,2);

accel_na = vec(ismrm_calculate_noise_amplification(unmixing, noise_matrix)); 
full_na  = vec(ismrm_calculate_noise_amplification(ccm, noise_matrix));

nonzero_ind = full_na > 0;

gmap = zeros(size(full_na));
gmap(nonzero_ind) = accel_na(nonzero_ind) ./ (full_na(nonzero_ind) .* acc_factor); 

gmap = reshape(gmap, [nx ny]);