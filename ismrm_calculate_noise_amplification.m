function na = ismrm_calculate_noise_amplification(unmixing, noise_matrix)
%
%  na = ismrm_calculate_noise_amplification(unmixing, noise_matrix)
%
%  Computes noise amplification from separate channel-by-channel images to
%  a combined single channel image.
%
%  INPUT:
%    - unmixing [x,y, coil]     : unmixing images for accelerated case,
%                                 of channel combination maps for the
%                                 unaccelerated case.
%    - noise_matrix [coil,coil] : noise covariance matrix
%
%  OUTPUT:
%    - na [x,y]                 : noise amplification map
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%
nx = size(unmixing, 1);
ny = size(unmixing, 2);
nc = size(unmixing, 3);

na = zeros([nx ny]);
for ic2 = 1:nc,
    for ic1 = 1:nc,
        na = na + noise_matrix(ic1, ic2) .* unmixing(:,:,ic1) .* conj(unmixing(:,:,ic2));
    end
end

na = sqrt(abs(na));