function unmix = ismrm_calculate_jer_unmixing(jer_lookup, acc_factor, ccm, regularization_scale, verbose)
%
%   unmix = ismrm_calculate_jer_unmixing(jer_lookup, acc_factor, ccm, regularization_scale, verbose)
%   
%   Calculates channel-by-channel local k-space unaliasing kernels based on
%   the provided joint-encoding relations and acceleration factor.
%
%   Transforms these kernels to image space and merges them with the
%   provided channel combination maps to create unmixing images.
%
%   INPUT:
%       jer_lookup [kx,ky,coil, kx, ky, coil] : Lookup table of joint
%                                               encoding relations. Kernel
%                                               extent take from jer size.
%       acc_factor  scalar          : Acceleration factor, e.g. 2
%       ccm         [x,y,coil]      : Channel combination maps
%       regularization_scale scalar : Controls aggressiveness of Tychonov 
%                                     regularization during calculation of 
%                                     unaliasing kernels.
%                                     0 = no regularization;
%                                     0.001 = default;
%                                     higher for more aggressive
%                                     regularization.
%       verbose     bool           : Set true for verbose output
%
%   OUTPUT:
%       unmix [x,y,coil]           : Image unmixing coefficients
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%


%%
% Compute unaliasing kernels
    
if (verbose),
    fprintf('Calculating unaliasing kernels...\n');
end
%tic
kernel_size = [size(jer_lookup,1) size(jer_lookup,2)];
target_location = bitshift(kernel_size, -1)+1;
num_channels = size(ccm, 3);
kernel = zeros(kernel_size(1),kernel_size(2), num_channels, num_channels);

for ic = 1:num_channels,
    kernel(target_location(1), target_location(2), ic, ic) = 1;
end

for s=1:(acc_factor-1),
   kernel_mask = zeros(kernel_size);
   kernel_mask(:,s:acc_factor:end) = 1;
   k = ismrm_compute_kspace_unaliasing_coefficients(jer_lookup, kernel_mask, regularization_scale);
   kernel = kernel + k;
end
%toc

%%
% Form unmixing images from channel combination maps and kernels
if (verbose),
    fprintf('Merging unaliasing and channel combination images...\n');
end
%tic
unmix = ismrm_calculate_unmixing_images_from_kspace_kernels(kernel, ccm);
%toc

if (verbose),
    fprintf('done.\n');
end

return
        
