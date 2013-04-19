function kernel = ismrm_compute_kspace_unaliasing_coefficients(jer_lookup, kernel_mask, regularization_scale)
%
%   kernel = ismrm_compute_kspace_unaliasing_coefficients(jer_lookup, kernel_mask, regularization_scale)
%   
%   Compute kspace unaliasing coefficients from joint encoding relations.
%
%
%   INPUT:
%       jer_lookup [kx,ky,coil, kx, ky, coil] : joint encoding relations lookup table
%       kernel_mask [kx,ky]        : e.g [1 1 1; 0 0 0; 1 1 1] for a 3x3
%                                    kernel with an acceleration factor of
%                                    2.
%       regularization_scale       : amount of Tychonov regularization to
%                                    apply.
%                                    0 (default) = no regularization.
%                                    0.001 moderate regularization
%                                    higher value for more regularization.
%
%
%   OUTPUT:
%       kernel [kx, ky, num_source_coils, num_target_coils] : k-space
%                 unaliasing kernels (for uniform undersampling pattern)
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%
if nargin < 3,
    regularization_scale = [];
end

if isempty(regularization_scale)
    regularization_scale = 0.0;
end

assert( regularization_scale >= 0, 'regularization_scale must be positive');

[kx_offsets, ky_offsets] = ind2sub(size(kernel_mask), find(kernel_mask == 1));

num_source = length(kx_offsets);
num_channel = size(jer_lookup,5);

Rss = zeros(num_source, num_channel, num_source, num_channel);
Rst = zeros(num_source, num_channel, num_channel);


for is2 = 1:num_source,
    for is1 = 1:num_source,
        Rss(is1, :, is2, :) = jer_lookup(kx_offsets(is1), ...
                                         ky_offsets(is1), ...
                                         kx_offsets(is2), ...
                                         ky_offsets(is2), ...
                                         :, :);
    end
end



kx_target = bitshift(size(kernel_mask,1),-1)+1;
ky_target = bitshift(size(kernel_mask,2),-1)+1;


for is1 = 1:num_source,
    Rst(is1, :, :) = jer_lookup(kx_offsets(is1), ...
                                ky_offsets(is1), ...
                                kx_target, ...
                                ky_target, ...
                                :, :);
            
end


num_basis = num_source * num_channel;
Rss = reshape(Rss, [num_basis num_basis]);

%%
%svals = svd(Rss);
%svals = svals / max(svals);
%regThreshold = ones(size(svals)) .* regularization_scale * mean(svals);
%figure; plot(svals);
%hold on; plot(regThreshold, 'r'); title('eigenvalues');

%figure; plot(1./ svals);
%hold on; plot(1./ (svals + regThreshold), 'r'); title('amplification');

%%
Rst = reshape(Rst, [num_basis num_channel]);

w = (Rss + eye(num_basis).* (regularization_scale * trace(Rss) / num_basis) ) \ Rst;

kernel = repmat(kernel_mask,[1 1 num_channel num_channel]);
kernel(kernel == 1) = w(:);
