function [unmix, gmap] = ismrm_calculate_sense_unmixing(acc_factor, csm, noise_matrix, regularization_factor)
%
% function [unmix, gmap] = ismrm_calculate_sense_unmixing(acc_factor, csm, noise_matrix, regularization_factor)
%
% Calculates the unmixing coefficients for a 2D image
%
% INPUT:
%       acc_factor  scalar             : Acceleration factor, e.g. 2
%       csm         [x, y, coil]       : Coil sensitivity map 
%       noise_matrix [coil, coil]      : noise covariance matrix
%       regularization_factor scaler   : adds Tychonov regularization.
%                                        0 = no regularization
%                                        0.001 = default
%                                        set higher for more aggressive
%                                        regularization.
%
% OUTPUT:
%       unmix       [x, y, coil] : Image unmixing coefficients for a single x location 
%       gmap        [x, y]       : Noise enhancement map 
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip Beatty (philip.beatty@sri.utoronto.ca)
%

assert(nargin > 1, 'At least 2 arguments needed');
assert(ndims(csm)==3, 'coil sensitivity map must have 3 dimensions');

if nargin < 3,
    noise_matrix = [];
end

if nargin < 4,
    regularization_factor = [];
end

if isempty(noise_matrix),
    noise_matrix = eye(size(csm,3));
end
if isempty(regularization_factor),
    regularization_factor = 0.00;
end

unmix = zeros(size(csm));

noise_matrix_inv = pinv(noise_matrix);

for x=1:size(csm,1), 
    unmix(x,:,:) = ismrm_calculate_sense_unmixing_1d(acc_factor, squeeze(csm(x,:,:)), noise_matrix_inv, regularization_factor); 
end

if (nargout > 1),
   gmap = sqrt(sum(abs(unmix).^2,3)).*sqrt(sum(abs(csm).^2,3));
end

function [unmix1d] = ismrm_calculate_sense_unmixing_1d(acc_factor, csm1d, noise_matrix_inv, regularization_factor)

[ny, nc] = size(csm1d);

if mod(ny, acc_factor) ~= 0,
    error('ny must be a multiple of acc_factor');
end

unmix1d = zeros(ny, nc);

n_blocks = ny/acc_factor;
for index = 1:n_blocks
    A = csm1d(index:n_blocks:ny,:).';
    if max(abs(A(:))) > 0,  
%        unmix1d(index:n_blocks:ny, :) = pinv(A);
        AHA = A'*noise_matrix_inv * A;
        reduced_eye = diag(abs(diag(AHA))>0);
        n_alias = sum(reduced_eye(:));
        scaled_reg_factor = regularization_factor * trace(AHA)/n_alias;
        
        unmix1d(index:n_blocks:ny, :) = pinv(AHA + reduced_eye .* scaled_reg_factor) * A' * noise_matrix_inv;
    end
end
